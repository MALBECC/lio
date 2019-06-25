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
!   * cdft_add_energy: adds CDFT terms to energy.                              !
!                                                                              !
! The following subroutines are only called internally:                        !
!   * cdft_get_chrg_constraint: Gets sum(atomic_charges) - charge_target.      !
!   * cdft_get_spin_constraint: Gets sum(atomic_spins) - spin_target.          !
!                                                                              !
!------------------------------------------------------------------------------!
! How to input CDFT options:                                                   !
!   CDFT options must be provided in the LIO input file (lio.in) but outside   !
!   namelists. Said file should contain the following sections:                !
!                                                                              !
!{CDFT}                                                                        !
!  CONST_CHARGE CONST_SPIN N_REG                                               !
!  REGION1_NATOM REGION1_CHARGE REGION1_SPIN                                   !
!  REGION2_NATOM REGION2_CHARGE REGION2_SPIN                                   !
!  REGIONN_NATOM REGIONN_CHARGE REGIONN_SPIN                                   !
!  REGION1_ATOMS                                                               !
!  REGION2_ATOMS                                                               !
!  REGIONN_ATOMS                                                               !
!{END}                                                                         !
!                                                                              !
! The {CDFT} and {END} terms indicate the beginning and end of CDFT input.     !
! N_REG (integer) is the number of regions to constrain, while CONST_CHARGE    !
! and CONST_SPIN (both integers too) indicate whether either/both charge       !
! (CONST_CHARGE=1) or/and spin (CONST_SPIN=1) are constrained.                 !
! After that first line, the following lines contain region information; each  !
! line belongs to a specific region of the molecule. REGION_CHARGE indicates   !
! the target charge for a region (in double precision), REGION_SPIN indicates  !
! the target spin of a region, and REGION_NATOM indicates the number of atoms  !
! in a region. Finally the last lines contain the atom indexes belonging to a  !
! region, in the same order as specified in the above lines. These should be   !
! written as a space-separated integer list. See the CDFT test in /tests for   !
! an example.                                                                  !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module cdft_data
   implicit none
   logical      :: cdft_chrg       = .false.
   logical      :: cdft_spin       = .false.


   ! Arrayable
   integer      :: cdft_atoms(200)
   real(kind=8) :: cdft_chrg_value = 0.1D0
   real(kind=8) :: cdft_spin_value = 0.0D0

   ! Other external vars.
   logical      :: doing_cdft = .true.

   ! Internal vars.
   real(kind=8) :: cdft_Vc     = 0.0D0  ! Potential for charge contraints.
   real(kind=8) :: cdft_Vs     = 0.0D0  ! Potential for spin constraints.
   real(kind=8) :: delta_V     = 0.0D0  ! Used for derivatives perturbation.
   real(kind=8) :: cdft_Vc_old = 0.0D0  ! For propagation.
   real(kind=8) :: cdft_Vs_old = 0.0D0  ! For propagation.
   real(kind=8) :: cdft_Cc_old = 0.0D0  ! For propagation.
   real(kind=8) :: cdft_Cs_old = 0.0D0  ! For propagation.
   real(kind=8) :: jacob(2,2)  = 0.0D0  ! Jacobian elements for advancing Vi
   type cdft_list
      integer      :: n_regions    = 0

      ! Main data for regions.
      integer     , allocatable :: r_natom(:)   ! Number of atoms.
      integer     , allocatable :: r_atoms(:,:) ! Atom indexes.
      real(kind=8), allocatable :: r_chrg(:)    ! Charge targets.
      real(kind=8), allocatable :: r_spin(:)    ! Spin targets.
      
      ! Atomic charges and spins.
      real(kind=8), allocatable :: at_chrg(:)
      real(kind=8), allocatable :: at_spin(:)
   end type cdft_list
   type(cdft_list) :: cdft_lists
end module cdft_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module cdft_subs
   implicit none
contains

! Checks if doing CDFT and sets energy_all_iterations to true in order to
! also calculate Becke charges in each iteration step.
subroutine cdft_options_check(do_becke)
   use cdft_data, only: cdft_chrg, cdft_spin, doing_cdft

   implicit none
   logical, intent(inout) :: do_becke

   if ((cdft_chrg) .or. (cdft_spin)) then
      doing_cdft     = .true.
      do_becke       = .true.
   endif
end subroutine cdft_options_check

! Reads a CDFT input from input file.
subroutine cdft_input_read(input_UID)
   use cdft_data, only: cdft_spin, cdft_chrg, cdft_lists
   implicit none
   integer, intent(in) :: input_UID
   
   character(len=10) :: buffer
   integer           :: ios, inp_spin, inp_chrg, ii
   
   rewind(input_UID)
   ios = 0
   do while ((trim(buffer) /= "{CDFT}") .and. (ios == 0) )
      read(input_UID,'(A10)', iostat=ios) buffer
   enddo

   ! If ios < 0, found EOF. No CDFT input provided.
   if (ios < 0) return

   ! Starts reading CDFT data.
   read(input_UID,*) inp_chrg, inp_spin, cdft_lists%n_regions
   if (inp_chrg == 1) cdft_chrg = .true.
   if (inp_spin == 1) cdft_spin = .true.

   if (allocated(cdft_lists%r_chrg))  deallocate(cdft_lists%r_chrg)
   if (allocated(cdft_lists%r_spin))  deallocate(cdft_lists%r_spin)
   if (allocated(cdft_lists%r_natom)) deallocate(cdft_lists%r_natom)
   allocate(cdft_lists%r_chrg(cdft_lists%n_regions))
   allocate(cdft_lists%r_spin(cdft_lists%n_regions))
   allocate(cdft_lists%r_natom(cdft_lists%n_regions))

   do ii = 1, cdft_lists%n_regions
      read(input_UID,*) cdft_lists%r_natom(ii), cdft_lists%r_chrg(ii), &
                        cdft_lists%r_spin(ii)
   enddo
   if (allocated(cdft_lists%r_atoms)) deallocate(cdft_lists%r_atoms)
   allocate(cdft_lists%r_atoms(cdft_lists%n_regions, &
                               maxval(cdft_lists%r_natom,1)))
   cdft_lists%r_atoms = 0
   do ii = 1, cdft_lists%n_regions
      read(input_UID,*) cdft_lists%r_atoms(ii,1:cdft_lists%r_natom(ii))
   enddo
end subroutine cdft_input_read

! Initializes arrays.
subroutine cdft_initialise(n_atoms)
   use cdft_data, only: cdft_lists, doing_cdft, cdft_atoms, cdft_Vc
   implicit none
   integer, intent(in) :: n_atoms

   if (allocated(cdft_lists%at_chrg)) deallocate(cdft_lists%at_chrg)
   if (allocated(cdft_lists%at_spin)) deallocate(cdft_lists%at_spin)
   allocate(cdft_lists%at_chrg(n_atoms))
   allocate(cdft_lists%at_spin(n_atoms))

   call g2g_cdft_init(doing_cdft, cdft_Vc, cdft_atoms)
end subroutine cdft_initialise

! Deallocates arrays.
subroutine cdft_finalise()
   use cdft_data, only: cdft_lists
   implicit none
   if (allocated(cdft_lists%at_chrg)) deallocate(cdft_lists%at_chrg)
   if (allocated(cdft_lists%at_spin)) deallocate(cdft_lists%at_spin)
end subroutine cdft_finalise

! Performs the CDFT iterative procedure.
! First, it gets dW/dVk (W being total energy and Vk the constrained
! potential) which equals the constraint value. Then stores the derivative
! and performs another cycle, extrapolating a new Vk from previous iterations.
subroutine cdft_set_potential(n_iter)
   use cdft_data, only: cdft_chrg, cdft_spin, cdft_Vc, cdft_Vs,      &
                        cdft_Vc_old, cdft_Vs_old, jacob, &
                        cdft_Cc_old, cdft_Cs_old
   implicit none
   integer     , intent(inout) :: n_iter

   real(kind=8) :: jacob_const, inv_jacob(2,2)
    
   ! Progates the potentials Vc/Vs using the inverse Jacobian.
   if (cdft_chrg .and. cdft_spin) then
      jacob_const = 1 / (jacob(2,2) * jacob(1,1) - jacob(2,1) * jacob(1,2))
      inv_jacob(1,1) =   jacob_const * jacob(2,2)
      inv_jacob(2,2) =   jacob_const * jacob(1,1)
      inv_jacob(2,1) = - jacob_const * jacob(1,2)
      inv_jacob(1,2) = - jacob_const * jacob(2,1)
   endif

   if (cdft_chrg) then
      if (cdft_spin) then
         cdft_Vc = cdft_Vc_old - 0.5D0 * (cdft_Cc_old * inv_jacob(1,1) + &
                                          cdft_Cs_old * inv_jacob(1,2))
      else
         cdft_Vc = cdft_Vc_old - 0.5D0 * cdft_Cc_old / jacob(1,1)
      endif
   endif

   if (cdft_spin) then
      if (cdft_chrg) then
         cdft_Vs = cdft_Vs_old - 0.5D0 * (cdft_Cc_old * inv_jacob(2,1) + &
                                          cdft_Cs_old * inv_jacob(2,2))
      else
         cdft_Vs = cdft_Vs_old - 0.5D0 * cdft_Cs_old / jacob(2,2)
      endif
      
   endif

   print*, "CDFT Vc: ", cdft_vc, "CDFT Vs: ", cdft_vs
end subroutine cdft_set_potential

! Makes a small perturbation in either Vc, Vs or both.
! This is used to approximate the Jacobian elements for
! Vi advancement in the iteration.
subroutine cdft_perturbation(mode)
   use cdft_data, only: cdft_Vs, cdft_Vc, jacob, delta_V, &
                        cdft_Vs_old, cdft_Vc_old
   implicit none
   integer, intent(in) :: mode

   ! Charge perturbation.
   if (mode == 1) then
      jacob(1,1)  = cdft_get_chrg_constraint()
      jacob(2,1)  = cdft_get_spin_constraint()
      cdft_Vs = cdft_Vs_old
      if (abs(cdft_Vc) < 1.0D-6) then
         delta_V = 1.0D-6
      else
         delta_V = abs(cdft_Vc * 0.01D0)
      endif
      cdft_Vc = cdft_Vc + delta_V
   endif

   ! Spin perturbation.
   if (mode == 2) then      
      jacob(1,2) = cdft_get_chrg_constraint()
      jacob(2,2) = cdft_get_spin_constraint()
      cdft_Vc = cdft_Vc_old
      if (abs(cdft_Vs) < 1.0D-6) then
         delta_V = 1.0D-6
      else
         delta_V = abs(cdft_Vs * 0.01D0)
      endif
      cdft_Vs = cdft_Vs + delta_V
   endif
end subroutine cdft_perturbation

subroutine cdft_set_jacobian(mode)
   use cdft_data, only: cdft_Vs, cdft_Vc, delta_V, jacob, &
                        cdft_Vs_old, cdft_Vc_old
   implicit none
   integer, intent(in) :: mode

   if (mode == 1) then
      jacob(1,1) = (cdft_get_chrg_constraint() - jacob(1,1)) / delta_V
      jacob(2,1) = (cdft_get_spin_constraint() - jacob(2,1)) / delta_V
      cdft_Vc    = cdft_Vc_old
   endif
   if (mode == 2) then
      jacob(1,2) = (cdft_get_chrg_constraint() - jacob(1,2)) / delta_V
      jacob(2,2) = (cdft_get_spin_constraint() - jacob(2,2)) / delta_V
      cdft_Vs    = cdft_Vs_old
   endif
end subroutine cdft_set_jacobian

! Checks if CDFT converged.
subroutine cdft_check_conver(rho_new, rho_old, converged)
   use cdft_data, only: cdft_lists, cdft_Cc_old, cdft_Cs_old, &
                        cdft_Vc_old, cdft_Vs_old, cdft_Vc, cdft_Vs
   implicit none
   real(kind=8), intent(in)  :: rho_new(:), rho_old(:)
   logical     , intent(out) :: converged 

   real(kind=8) :: rho_diff
   integer      :: jj
   
   rho_diff = 0.0D0
   do jj = 1 , size(rho_new,1)
      rho_diff  = rho_diff + (rho_new(jj) - rho_old(jj)) * &
                             (rho_new(jj) - rho_old(jj))
   enddo
   rho_diff = sqrt(rho_diff) / dble(size(rho_new,1))
   cdft_Vc_old = cdft_Vc
   cdft_Vs_old = cdft_Vs
   cdft_Cc_old = cdft_get_chrg_constraint()
   cdft_Cs_old = cdft_get_spin_constraint()

   write(*,*) "CDFT convergence:"
   write(*,*) "Rho: ", rho_diff, "Constraint: ", cdft_Cc_old, cdft_Cs_old
   converged = .false.
   if ((rho_diff < 1D-4) .and. (abs(cdft_Cc_old) < 1D-9)) converged = .true.
end subroutine cdft_check_conver

! Adds CDFT terms to total energy.
subroutine cdft_add_energy(energ, mode_in)
   ! MODES:
   !   0 = add both spin and charge constraint energy (if able).
   !   1 = add only charge.
   !   2 = add only spin.
   use cdft_data, only: cdft_chrg, cdft_spin, cdft_Vc, cdft_Vs
   implicit none
   integer     , intent(in), optional :: mode_in
   real(kind=8), intent(inout)        :: energ

   real(kind=8) :: c_value = 0.0D0
   integer      :: mode    = 0
   
   if (present(mode_in)) mode = mode_in

   ! First adds charge constraints.
   if (((cdft_chrg) .and. (mode == 0)) .or. (mode == 1)) then
      c_value = cdft_get_chrg_constraint()
      energ   = energ + c_value * cdft_Vc
   endif

   ! Then adds spin constraints.
   if (((cdft_spin) .and. (mode == 0)) .or. (mode == 2)) then
      c_value = cdft_get_spin_constraint()
      energ   = energ + c_value * cdft_Vs
   endif
end subroutine

! Gets Becke atomic charges and calculates its difference
! with the constrained total value.
function cdft_get_chrg_constraint() result(const_value)
   use cdft_data, only: cdft_atoms, cdft_chrg_value, cdft_lists
   implicit none
   integer                   :: c_index
   real(kind=8)              :: const_value

   call g2g_get_becke_dens(cdft_lists%at_chrg)

   c_index     = 1
   const_value = 0.0D0
   do while (cdft_atoms(c_index) > 0)
      const_value = const_value + cdft_lists%at_chrg(cdft_atoms(c_index))
      c_index = c_index +1
   enddo
   const_value = cdft_chrg_value - const_value

   return
end function cdft_get_chrg_constraint

! Gets Becke atomic spin and calculates its difference
! with the constrained total value.
function cdft_get_spin_constraint() result(const_value)
   use cdft_data, only: cdft_atoms, cdft_spin_value, cdft_lists
   implicit none
   integer                    :: c_index
   real(kind=8)               :: const_value

   call g2g_get_becke_spin(cdft_lists%at_spin)

   c_index     = 1
   const_value = 0.0D0
   do while (cdft_atoms(c_index) > 0)
      const_value = const_value + cdft_lists%at_spin(cdft_atoms(c_index))
      c_index = c_index +1
   enddo
   const_value = cdft_spin_value - const_value

   return
end function cdft_get_spin_constraint

end module cdft_subs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
