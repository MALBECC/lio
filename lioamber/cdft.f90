!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% CONSTRAINED DFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains routines (cdft_subs) and variables (cdft_data) for        !
! Constrained DFT (CDFT) calculations. This implementation uses the Becke      !
! partitioning scheme for the constraints, as used by Holmberg and Laasoen     !
! (JCTC 2016, DOI: 10.1021/acs.jctc.6b01085).                                  !
!                                                                              !
! Important input variables are:                                               !
!   * cdft_atoms: integer array, indicates which atoms are participatin in the !
!                 constrain.                                                   !
!   * cdft_chrg: logical, indicates if atomic charge is constrained.           !
!   * cdft_spin: logical, indicates if atomic spin is constrained.             !
!   * cdft_chrg_value: double precision real, indicates the value for the atom !
!                      charge constrain.                                       !
!   * cdft_spin_value: double precision real, indicates the value for the atom !
!                      spin constrain.                                         !
!                                                                              !
! The following subroutines are called externally:                             !
!   * cdft_check_options: checks internal consistency in input options.        !
!   * cdft_add_fock: adds CDFT terms to Fock matrix.                           !
!   * cdft_add_energy: adds CDFT terms to energy.                              !
!                                                                              !
! The following subroutines are only called internally:                        !
!   * cdft_get_chrg_constraint: Gets sum(atomic_charges) - charge_target.      !
!   * cdft_get_spin_constraint: Gets sum(atomic_spins) - spin_target.          !
!                                                                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module cdft_data
   implicit none
   logical      :: cdft_chrg       = .false.
   logical      :: cdft_spin       = .false.
   integer      :: cdft_atoms(200) = 0
   real(kind=8) :: cdft_chrg_value = 0.0D0
   real(kind=8) :: cdft_spin_value = 0.0D0

   ! Other external vars.
   logical      :: doing_cdft = .false.

   ! Internal vars.
   real(kind=8) :: cdft_cV = 1.0D0  ! Potential for charge contraints.
   real(kind=8) :: cdft_sV = 1.0D0  ! Potential for spin constraints.
   integer      :: n_cyc   = 0
   real(kind=8) :: cdft_cV_list(10) = 0.0D0
   real(kind=8) :: cdft_sV_list(10) = 0.0D0
end module cdft_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module cdft_subs
   implicit none
contains

! Performs the CDFT iterative procedure.
subroutine cdft_set_potential(energ, atom_z, n_atoms)
   use cdft_data, only: cdft_chrg, cdft_spin
   implicit none
   integer     , intent(in)    :: n_atoms, atom_z(n_atoms)
   real(kind=8), intent(inout) :: energ

   real(kind=8) :: c_value, energ_tmp

   if (cdft_chrg) then
      energ_tmp = 0.0D0
      call cdft_add_energy(energ_tmp, atom_z, n_atoms, 1)
      
   endif

   if (cdft_spin) then
      energ_tmp = 0.0D0
      call cdft_add_energy(energ_tmp, atom_z, n_atoms, 2)
      
   endif
end subroutine cdft_set_potential


! Checks if CDFT converged.
subroutine cdft_check_conver(rho_new, rho_old, converged)
   implicit none
   real(kind=8), intent(in)  :: rho_new(:), rho_old(:)
   logical     , intent(out) :: converged 

   real(kind=8) :: rho_diff
   integer      :: jj
   
   rho_diff = 0.0D0
   do jj = 1 , size(rho_new,1)
      rho_diff  = rho_diff + rho_new(jj) - rho_old(jj) * &
                             rho_new(jj) - rho_old(jj)
   enddo
   rho_diff = sqrt(rho_diff) / dble(size(rho_new,1))

   converged = .false.
   if (rho_diff < 1D-5) converged = .true.

   do jj = 1, 9
      cdft_cV_list(jj) = cdft_cV_list(jj+1)
      cdft_sV_list(jj) = cdft_sV_list(jj+1)
   enddo
   cdft_cV_list(10) == cdft_cV
   cdft_sV_list(10) == cdft_sV

end subroutine cdft_check_conver

! Checks if doing CDFT and sets energy_all_iterations to true in order to
! also calculate Becke charges in each iteration step.
subroutine cdft_check_options(energ_all_iter)
   use cdft_data, only: cdft_chrg, cdft_spin, doing_cdft

   implicit none
   logical, intent(inout) :: energ_all_iter

   if ((cdft_chrg) .or. (cdft_spin)) then
      doing_cdft     = .true.
      energ_all_iter = .true.
   endif
end subroutine cdft_check_options

! Adds CDFT potential terms to Fock matrix in atomic orbital basis.
subroutine cdft_add_fock(Fmat, Smat, atom_z, n_atoms)
   use cdft_data, only: cdft_cV, cdft_sV, cdft_spin, cdft_chrg
   implicit none
   integer     , intent(in)    :: n_atoms, atom_z(n_atoms)
   real(kind=8), intent(in)    :: Smat(:,:)
   real(kind=8), intent(inout) :: Fmat(:,:)

   real(kind=8) :: c_value = 0.0D0

   ! First adds charge constraints.
   if (cdft_chrg) then
      c_value = cdft_get_chrg_constraint(atom_z, n_atoms)
      Fmat    = Fmat + Smat * c_value * cdft_cV
   endif

   ! Then adds spin constraints.
   if (cdft_spin) then
      c_value = cdft_get_spin_constraint(n_atoms)
      Fmat    = Fmat + Smat * c_value * cdft_sV
   endif

end subroutine cdft_add_fock

! Adds CDFT terms to total energy.
subroutine cdft_add_energy(energ, atom_z, n_atoms, mode_in)
   ! MODES:
   !   0 = add both spin and charge constraint energy (if able).
   !   1 = add only charge.
   !   2 = add only spin.
   implicit none
   integer     , intent(in), optional :: mode_in
   integer     , intent(in)           :: n_atoms, atom_z(n_atoms)
   real(kind=8), intent(inout)        :: energ

   real(kind=8) :: c_value = 0.0D0
   integer      :: mode    = 0
   
   if (present(mode_in)) mode = mode_in

   ! First adds charge constraints.
   if (((cdft_chrg) .and. (mode == 0)) .or. (mode == 1)) then
      c_value = cdft_get_chrg_constraint(atom_z, n_atoms)
      energ   = energ + c_value * cdft_cV
   endif

   ! Then adds spin constraints.
   if (((cdft_spin) .and. (mode == 0)) .or. (mode == 2)) then
      c_value = cdft_get_spin_constraint(n_atoms)
      energ   = energ + c_value * cdft_sV
   endif
end subroutine

! Gets Becke atomic charges and calculates its difference
! with the constrained total value.
function cdft_get_chrg_constraint(atom_z, n_atoms) result(const_value)
   use cdft_data, only: cdft_atoms, cdft_chrg_value
   implicit none
   integer, intent(in) :: n_atoms, atom_z(n_atoms)

   integer      :: c_index     = 1
   real(kind=8) :: const_value = 0.0D0
   real(kind=8) :: chrg(:)

   allocate(chrg(n_atoms))
   call g2g_cdft_get_becke_dens(chrg)
   chrg = chrg - dble(atom_z)

   do while (cdft_atoms(c_index) > 0)
      const_value = const_value + chrg(cdft_atoms(c_index))
      c_index = c_index +1
   enddo
   const_value = const_value - cdft_chrg_value

   deallocate(chrg)
   return
end function cdft_get_chrg_constraint

! Gets Becke atomic spin and calculates its difference
! with the constrained total value.
function cdft_get_spin_constraint(n_atoms) result(const_value)
   use cdft_data, only: cdft_atoms, cdft_spin_value
   implicit none
   integer, intent(in) :: n_atoms

   integer      :: c_index     = 1
   real(kind=8) :: const_value = 0.0D0
   real(kind=8) :: spin(:)

   allocate(spin(n_atoms))
   call g2g_cdft_get_becke_dens(spin)

   do while (cdft_atoms(c_index) > 0)
      const_value = const_value + spin(cdft_atoms(c_index))
      c_index = c_index +1
   enddo
   const_value = const_value - cdft_spin_value

   deallocate(spin)
   return
end function cdft_get_spin_constraint


end module cdft_subs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!