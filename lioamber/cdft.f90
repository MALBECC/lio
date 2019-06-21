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
   logical      :: cdft_chrg       = .true.
   logical      :: cdft_spin       = .false.
   integer      :: cdft_atoms(200)
   real(kind=8) :: cdft_chrg_value = 1.0D0
   real(kind=8) :: cdft_spin_value = 0.0D0

   ! Other external vars.
   logical      :: doing_cdft = .true.
   real(kind=8) :: HL_gap     = 0.0D0

   ! Internal vars.
   real(kind=8) :: cdft_Vc = 0.0D0  ! Potential for charge contraints.
   real(kind=8) :: cdft_Vs = 0.0D0  ! Potential for spin constraints.

   type cdft_list
      integer      :: list_size = 10
      real(kind=8) :: Vc(10)    = 0.0D0 ! Potential for charge contraints.
      real(kind=8) :: Vs(10)    = 0.0D0 ! Potential for spin contraints.
      real(kind=8) :: dWdVc(10) = 0.0D0 ! First energy derivative from Vc.
      real(kind=8) :: dWdVs(10) = 0.0D0 ! First energy derivative from Vs.
   end type cdft_list
   type(cdft_list) :: cdft_lists
end module cdft_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module cdft_subs
   implicit none
contains

! Performs the CDFT iterative procedure.
! First, it gets dW/dVk (W being total energy and Vk the constrained
! potential) which equals the constraint value. Then stores the derivative
! and performs another cycle, extrapolating a new Vk from previous iterations.
subroutine cdft_set_potential(n_iter, atom_z, n_atoms)
   use cdft_data, only: cdft_chrg, cdft_spin, cdft_Vc, cdft_Vs, HL_gap, &
                        cdft_lists
   implicit none
   integer     , intent(in)    :: n_atoms, atom_z(n_atoms)
   integer     , intent(inout) :: n_iter

   integer      :: icount, start_ind
   real(kind=8) :: c_value, first_dev

   if (n_iter == 1) then
      cdft_Vc = 0.0D0
      cdft_Vs = 0.0D0
      return
   else if (n_iter == 2) then
      cdft_Vc = 1.0D10
      cdft_Vs = 1.0D10
   else 
      ! Third iteration
      cdft_Vc = cdft_Vc * -1.0D0
      cdft_Vs = cdft_Vs * -1.0D0
   endif
      
   ! Calculates the first derivatives.
   ! This start index is to avoid resizing or recalculating indexes in each 
   ! iteration. It starts from +2 since the current derivatives correspond
   ! to the previous iteration step (i.e., when n_iter=2 we only use one point
   ! since data corresponds only the first iteration).
   start_ind = max(1, cdft_lists%list_size - n_iter +2)
   if (cdft_chrg) then
      cdft_lists%dWdVc(cdft_lists%list_size) = &
                                      cdft_get_chrg_constraint(atom_z, n_atoms)
      if (n_iter > 3) then
         cdft_Vc = cdft_max_V(cdft_lists%Vc(start_ind:cdft_lists%list_size), &
                              cdft_lists%dWdVc(start_ind:cdft_lists%list_size))
      endif
   endif

   if (cdft_spin) then
      cdft_lists%dWdVs(cdft_lists%list_size) = &
                                      cdft_get_spin_constraint(n_atoms)
      if (n_iter > 3) then
         cdft_Vs = cdft_max_V(cdft_lists%Vs(start_ind:cdft_lists%list_size), &
                              cdft_lists%dWdVs(start_ind:cdft_lists%list_size))
      endif
   endif

   ! Moves the list of previous values for Vk.
   do icount = 1, cdft_lists%list_size -1
      cdft_lists%Vc(icount)    = cdft_lists%Vc(icount+1)
      cdft_lists%Vs(icount)    = cdft_lists%Vs(icount+1)
      cdft_lists%dWdVc(icount) = cdft_lists%dWdVc(icount+1)
      cdft_lists%dWdVs(icount) = cdft_lists%dWdVs(icount+1)
   enddo
   cdft_lists%Vc(cdft_lists%list_size) = cdft_Vc
   cdft_lists%Vs(cdft_lists%list_size) = cdft_Vs 

   print*, "CDFT Vc: ", cdft_vc, "CDFT Vs: ", cdft_vs
end subroutine cdft_set_potential

! Sets the HOMO-LUMO gap, which is really used only the first time
! as an estimator for the Vk used as constraints potentials.
subroutine cdft_set_HL_gap(HL_in)
   use cdft_data, only: HL_gap
   implicit none
   real(kind=8), intent(in) :: HL_in

   HL_gap = HL_in
end subroutine cdft_set_HL_gap

! Checks if CDFT converged.
subroutine cdft_check_conver(rho_new, rho_old, converged)
   use cdft_data, only: cdft_Vc, cdft_Vs, cdft_lists
   implicit none
   real(kind=8), intent(in)  :: rho_new(:), rho_old(:)
   logical     , intent(out) :: converged 

   real(kind=8) :: rho_diff, Vc_diff, Vs_diff
   integer      :: jj
   
   rho_diff = 0.0D0
   do jj = 1 , size(rho_new,1)
      rho_diff  = rho_diff + (rho_new(jj) - rho_old(jj)) * &
                             (rho_new(jj) - rho_old(jj))
   enddo
   rho_diff = sqrt(rho_diff) / dble(size(rho_new,1))
   Vc_diff  = abs(cdft_lists%Vc(cdft_lists%list_size) &
                - cdft_lists%Vc(cdft_lists%list_size-1))
   Vs_diff  = abs(cdft_lists%Vs(cdft_lists%list_size) &
            - cdft_lists%Vs(cdft_lists%list_size-1))


   write(*,*) "CDFT convergence:"
   write(*,*) "Rho:", rho_diff, "Vc: ", Vc_diff, "Vs: ", Vs_diff
   converged = .false.
   if ((rho_diff < 1D-5) .and. (Vc_diff < 1D-4) ) converged = .true.
end subroutine cdft_check_conver

! Checks if doing CDFT and sets energy_all_iterations to true in order to
! also calculate Becke charges in each iteration step.
subroutine cdft_options_check(energ_all_iter, do_becke)
   use cdft_data, only: cdft_chrg, cdft_spin, doing_cdft

   implicit none
   logical, intent(inout) :: energ_all_iter, do_becke

   if ((cdft_chrg) .or. (cdft_spin)) then
      doing_cdft     = .true.
      energ_all_iter = .true.
      do_becke       = .true.
   endif
end subroutine cdft_options_check

! Adds CDFT potential terms to core Fock matrix in atomic orbital basis.
subroutine cdft_add_fock(Fmat, Smat, atom_z, n_atoms)
   use cdft_data, only: cdft_Vc, cdft_Vs, cdft_spin, cdft_chrg
   implicit none
   integer     , intent(in)    :: n_atoms, atom_z(n_atoms)
   real(kind=8), intent(in)    :: Smat(:,:)
   real(kind=8), intent(inout) :: Fmat(:,:)

   real(kind=8) :: c_value = 0.0D0
   integer :: ii, jj

   ! First adds charge constraints.
   if (cdft_chrg) then
      c_value = cdft_get_chrg_constraint(atom_z, n_atoms)
      do ii=1, size(Fmat,1)
      do jj=1, size(Fmat,2)
         Fmat(ii,jj)    = Fmat(ii,jj) + Smat(ii,jj) * c_value * cdft_Vc
      enddo
      enddo
   endif

   ! Then adds spin constraints.
   if (cdft_spin) then
      c_value = cdft_get_spin_constraint(n_atoms)
      Fmat    = Fmat + Smat * c_value * cdft_Vs
   endif

end subroutine cdft_add_fock

! Adds CDFT terms to total energy.
subroutine cdft_add_energy(energ, atom_z, n_atoms, mode_in)
   ! MODES:
   !   0 = add both spin and charge constraint energy (if able).
   !   1 = add only charge.
   !   2 = add only spin.
   use cdft_data, only: cdft_chrg, cdft_spin, cdft_Vc, cdft_Vs
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
      energ   = energ + c_value * cdft_Vc
   endif

   ! Then adds spin constraints.
   if (((cdft_spin) .and. (mode == 0)) .or. (mode == 2)) then
      c_value = cdft_get_spin_constraint(n_atoms)
      energ   = energ + c_value * cdft_Vs
   endif
end subroutine

! Gets Becke atomic charges and calculates its difference
! with the constrained total value.
function cdft_get_chrg_constraint(atom_z, n_atoms) result(const_value)
   use cdft_data, only: cdft_atoms, cdft_chrg_value
   implicit none
   integer     , intent(in)  :: n_atoms, atom_z(n_atoms)

   integer                   :: c_index
   real(kind=8)              :: const_value
   real(kind=8), allocatable :: chrg(:)

   allocate(chrg(n_atoms))
   call g2g_get_becke_dens(chrg)
   chrg = dble(atom_z) - chrg

   c_index     = 1
   const_value = 0.0D0
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
   integer, intent(in)        :: n_atoms

   integer                    :: c_index
   real(kind=8)               :: const_value
   real(kind=8), allocatable  :: spin(:)

   allocate(spin(n_atoms))
   call g2g_get_becke_spin(spin)

   c_index     = 1
   const_value = 0.0D0
   do while (cdft_atoms(c_index) > 0)
      const_value = const_value + spin(cdft_atoms(c_index))
      c_index = c_index +1
   enddo
   const_value = const_value - cdft_spin_value

   deallocate(spin)
   return
end function cdft_get_spin_constraint

function cdft_max_V(V, foV) result(max_V)
   implicit none
   real(kind=8), intent(in) :: V(:), foV(:)

   integer      :: max_ind, point(3)
   real(kind=8) :: max_V, f_a, f_b

   max_ind = maxloc(foV,1)
   if (size(foV,1) == 1) then
      max_V = V(1)
      return
   else if (size(foV,1) == 1) then
      max_V = V(max_ind)
      return
   endif

   ! Performs a quadratic inter/extrapolation in order to get
   ! the value of V for which foV is a maximum.
   ! First, we get a, b and c for foV = a.V^2 + b.V + c by means of
   ! ancient black magic:
   !     a) y1 = a * x1^2 + b * x1 + c
   !     b) y2 - y1 = a * (x2^2 - x1^2) + b * (x2 - x1)
   !        => (y2 - y1) / (x2 - x1) = b + a * (x2 + x1) = J21
   !     c) J31 - J21 = a * (x3 - x2)
   !        => a = (J31 - J21) / (x3 - x2)
   ! Points 1, 2 and 3 are chosen according to maxloc result.
   ! We can do this since for CDFT it seems W(Vk) is always concave
   ! (or convex, whichever it is).
   if (max_ind == 1) then
      point(1) = 3
      point(2) = 1
      point(3) = 2
   else if (max_ind == size(foV,1)) then
      point(1) = size(foV,1) -2
      point(2) = size(foV,1)
      point(3) = size(foV,1) -1
   else
      point(1) = max_ind -1
      point(2) = max_ind
      point(3) = max_ind +1
   endif

   f_a = (  (fov(point(3)) - fov(point(2))) / (V(point(3)) - V(point(2)))  &
          - (fov(point(1)) - fov(point(2))) / (V(point(1)) - V(point(2)))) &
          / (V(point(3)) - V(point(1)))
   f_b = (fov(point(3)) - fov(point(2))) / (V(point(3)) - V(point(2))) &
         - f_a * (V(point(3)) + V(point(2)))

   ! Finally, we obtain the position of the maximum by simple derivation.
   max_V = - 0.5D0 * f_b / f_a
   print*, "FOV", fov, "V", V
   print*, "Max_V", max_V, f_b, f_a
   return 
end function cdft_max_V


end module cdft_subs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!