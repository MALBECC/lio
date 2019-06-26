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
   logical      :: cdft_chrg  = .false.
   logical      :: cdft_spin  = .false.
   logical      :: doing_cdft = .false.

   real(kind=8), allocatable :: at_chrg(:)     ! List of atomic charges.
   real(kind=8), allocatable :: at_spin(:)     ! List of atomic spin charges.
   real(kind=8), allocatable :: jacob(:,:)     ! Jacobian for advancing Vi
   real(kind=8), allocatable :: jacob_i(:,:)   ! Inverse of the Jacobian
   real(kind=8)              :: j_step = 0.5D0 ! Step for Jacobian advancement.
   real(kind=8)              :: c_prev = 100.0 ! For convergence.

   type cdft_region_data
      integer      :: n_regions = 0 ! Number of regions.
      integer      :: max_nat   = 0 ! Maximum number of atoms in a region.

      ! Main data for regions.
      integer     , allocatable :: natom(:)   ! Number of atoms.
      integer     , allocatable :: atoms(:,:) ! Atom indexes.
      real(kind=8), allocatable :: chrg(:)    ! Charge targets.
      real(kind=8), allocatable :: spin(:)    ! Spin targets.
      real(kind=8), allocatable :: Vc(:)      ! Charge potentials.
      real(kind=8), allocatable :: Vs(:)      ! Spin potentials.
      real(kind=8), allocatable :: Cc(:)      ! Charge constraint value.
      real(kind=8), allocatable :: Cs(:)      ! Spin constraint value.
      
      ! Arrays for Jacobian calculation.
      real(kind=8), allocatable :: Vc_old(:) ! Charge potentials.
      real(kind=8), allocatable :: Vs_old(:) ! Spin potentials.
      real(kind=8), allocatable :: Cc_old(:) ! Charge constraint value.
      real(kind=8), allocatable :: Cs_old(:) ! Spin constraint value.

   end type cdft_region_data

   type(cdft_region_data) :: cdft_reg
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
   use cdft_data, only: cdft_spin, cdft_chrg, cdft_reg
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
   read(input_UID,*) inp_chrg, inp_spin, cdft_reg%n_regions
   if (inp_chrg == 1) cdft_chrg = .true.
   if (inp_spin == 1) cdft_spin = .true.

   if (allocated(cdft_reg%chrg))  deallocate(cdft_reg%chrg)
   if (allocated(cdft_reg%spin))  deallocate(cdft_reg%spin)
   if (allocated(cdft_reg%natom)) deallocate(cdft_reg%natom)
   allocate(cdft_reg%chrg(cdft_reg%n_regions))
   allocate(cdft_reg%spin(cdft_reg%n_regions))
   allocate(cdft_reg%natom(cdft_reg%n_regions))

   do ii = 1, cdft_reg%n_regions
      read(input_UID,*) cdft_reg%natom(ii), cdft_reg%chrg(ii), &
                        cdft_reg%spin(ii)
   enddo
   if (allocated(cdft_reg%atoms)) deallocate(cdft_reg%atoms)
   cdft_reg%max_nat = maxval(cdft_reg%natom,1)
   allocate(cdft_reg%atoms(cdft_reg%n_regions, cdft_reg%max_nat))
   cdft_reg%atoms = 0
   do ii = 1, cdft_reg%n_regions
      read(input_UID,*) cdft_reg%atoms(ii,1:cdft_reg%natom(ii))
   enddo
end subroutine cdft_input_read

! Initializes arrays.
subroutine cdft_initialise(n_atoms)
   use cdft_data, only: cdft_reg, cdft_chrg, cdft_spin, at_chrg, at_spin, &
                        jacob, jacob_i
   implicit none
   integer, intent(in) :: n_atoms
   integer             :: J_size = 0

   if (cdft_chrg .or. cdft_spin) then
      if (allocated(at_chrg))         deallocate(at_chrg)
      if (allocated(cdft_reg%Vc))     deallocate(cdft_reg%Vc)
      if (allocated(cdft_reg%Cc))     deallocate(cdft_reg%Cc)
      if (allocated(cdft_reg%Vc_old)) deallocate(cdft_reg%Vc_old)
      if (allocated(cdft_reg%Cc_old)) deallocate(cdft_reg%Cc_old)
      allocate(at_chrg(n_atoms))
      allocate(cdft_reg%Vc(cdft_reg%n_regions))
      allocate(cdft_reg%Cc(cdft_reg%n_regions))
      allocate(cdft_reg%Vc_old(cdft_reg%n_regions))
      allocate(cdft_reg%Cc_old(cdft_reg%n_regions))

      if (allocated(at_spin))         deallocate(at_spin)
      if (allocated(cdft_reg%Vs))     deallocate(cdft_reg%Vs)
      if (allocated(cdft_reg%Cs))     deallocate(cdft_reg%Cs)
      if (allocated(cdft_reg%Vs_old)) deallocate(cdft_reg%Vs_old)
      if (allocated(cdft_reg%Cs_old)) deallocate(cdft_reg%Cs_old)
      allocate(at_spin(n_atoms))  
      allocate(cdft_reg%Vs(cdft_reg%n_regions))
      allocate(cdft_reg%Cs(cdft_reg%n_regions))
      allocate(cdft_reg%Vs_old(cdft_reg%n_regions))
      allocate(cdft_reg%Cs_old(cdft_reg%n_regions))
   endif

   if (allocated(jacob)) deallocate(jacob)
   if (cdft_chrg .and. cdft_spin) then
      J_size = 2 * cdft_reg%n_regions
   else if (cdft_chrg .or. cdft_spin) then
      J_size = cdft_reg%n_regions
   endif
   allocate(jacob(J_size, J_size))

   if (J_size > 1) then
      if (allocated(jacob_i)) deallocate(jacob_i)
      allocate(jacob_i(J_size, J_size))
   endif

   call g2g_cdft_init(cdft_chrg, cdft_spin, cdft_reg%n_regions, &
                      cdft_reg%max_nat, cdft_reg%natom,     &
                      cdft_reg%atoms)
end subroutine cdft_initialise

! Deallocates arrays.
subroutine cdft_finalise()
   use cdft_data, only: cdft_reg, at_spin, at_chrg, jacob
   implicit none
   if (allocated(jacob))           deallocate(jacob)
   if (allocated(at_chrg))         deallocate(at_chrg)
   if (allocated(at_spin))         deallocate(at_spin)
   if (allocated(cdft_reg%Vc))     deallocate(cdft_reg%Vc)
   if (allocated(cdft_reg%Vs))     deallocate(cdft_reg%Vs)
   if (allocated(cdft_reg%natom))  deallocate(cdft_reg%natom)
   if (allocated(cdft_reg%atoms))  deallocate(cdft_reg%atoms)
   if (allocated(cdft_reg%chrg))   deallocate(cdft_reg%chrg)
   if (allocated(cdft_reg%spin))   deallocate(cdft_reg%spin)
   if (allocated(cdft_reg%Cc))     deallocate(cdft_reg%Cc)
   if (allocated(cdft_reg%Vc_old)) deallocate(cdft_reg%Vc_old)
   if (allocated(cdft_reg%Cc_old)) deallocate(cdft_reg%Cc_old)
   if (allocated(cdft_reg%Cs))     deallocate(cdft_reg%Cs)
   if (allocated(cdft_reg%Vs_old)) deallocate(cdft_reg%Vs_old)
   if (allocated(cdft_reg%Cs_old)) deallocate(cdft_reg%Cs_old)
end subroutine cdft_finalise

! Performs the CDFT iterative procedure.
subroutine CDFT(fock_a, rho_a, fock_b, rho_b, Pmat_vec, natom)
   use typedef_operator, only: operator

   implicit none
   integer       , intent(in)    :: natom
   real(kind=8)  , intent(inout) :: Pmat_vec(:)
   type(operator), intent(inout) :: fock_a, rho_a, fock_b, rho_b

   integer      :: cdft_iter, max_cdft_iter
   logical      :: cdft_converged = .false.
   real(kind=8) :: energ
   real(kind=8)  , allocatable :: Pmat_old(:)

   max_cdft_iter = 100
   cdft_iter     = 0
   allocate(Pmat_old(size(Pmat_vec,1)))
   call cdft_initialise(natom)
   do while ((.not. cdft_converged) .and. (cdft_iter < max_cdft_iter))
      cdft_iter = cdft_iter +1
      Pmat_old = Pmat_vec
      call SCF(energ, fock_a, rho_a, fock_b, rho_b)
      call cdft_check_conver(Pmat_vec, Pmat_old, cdft_converged, cdft_iter)

      if (.not. cdft_converged) then
         ! Calculates perturbations and Jacobian.
         call cdft_get_jacobian(fock_a, rho_a, fock_b, rho_b)
         call cdft_set_potential()
      endif
   enddo

   call cdft_finalise()
   deallocate(Pmat_old)
end subroutine CDFT

! Gets the jacobian matrix by making a small perturbation in each direction.
! Also, if more than one constraint is present, calculates the inverse
! jacobian by means of LAPACK's DGETRI.
subroutine cdft_get_jacobian(fock_a, rho_a, fock_b, rho_b)
   use typedef_operator, only: operator
   use cdft_data       , only: cdft_chrg, cdft_spin, cdft_reg, jacob, jacob_i
   
   implicit none
   type(operator), intent(inout) :: fock_a, rho_a, fock_b, rho_b
   integer      :: ii, jj
   real(kind=8) :: energ, dV

   ! Variables for LAPACK
   integer                   :: LWORK, INFO
   integer     , allocatable :: IPIV(:,:)
   real(kind=8), allocatable :: WORK(:)

   jacob = 0.0D0
   if (cdft_chrg) then
      cdft_reg%Vc_old = cdft_reg%Vc
      call cdft_get_chrg_constraints(cdft_reg%Cc_old)
      if (cdft_spin) call cdft_get_spin_constraints(cdft_reg%Cs_old)

      do ii = 1, cdft_reg%n_regions
         if (abs(cdft_reg%Vc(ii)) < 1.0D-8) then
            dV = 1.0D-9
         else
            dV = abs(cdft_reg%Vc(ii) * 0.1D0)
         endif
         cdft_reg%Vc(ii) = cdft_reg%Vc(ii) + dV

         call g2g_cdft_set_v(cdft_reg%Vc, cdft_reg%Vs)
         call SCF(energ, fock_a, rho_a, fock_b, rho_b)

         call cdft_get_chrg_constraints(cdft_reg%Cc)
         do jj = 1, cdft_reg%n_regions
            jacob(jj,ii) = jacob(jj,ii) + (cdft_reg%Cc(jj) - cdft_reg%Cc_old(jj)) / dV
         enddo
         if (cdft_spin) then
            call cdft_get_spin_constraints(cdft_reg%Cs)
            do jj = 1, cdft_reg%n_regions
               jacob(jj+cdft_reg%n_regions,ii) = &
                           (cdft_reg%Cs(jj) - cdft_reg%Cs_old(jj)) / dV
            enddo
         endif
         cdft_reg%Vc = cdft_reg%Vc_old   
      enddo
      cdft_reg%Cc = cdft_reg%Cc_old
      if (cdft_spin) cdft_reg%Cs = cdft_reg%Cs_old
   endif

   if (cdft_spin) then
      cdft_reg%Vs_old = cdft_reg%Vs
      call cdft_get_spin_constraints(cdft_reg%Cs_old)
      if (cdft_chrg) call cdft_get_chrg_constraints(cdft_reg%Cc_old)

      do ii = 1, cdft_reg%n_regions
         if (abs(cdft_reg%Vs(ii)) < 1.0D-8) then
            dV = 1.0D-10
         else
            dV = abs(cdft_reg%Vs(ii) * 0.01D0)
         endif
         cdft_reg%Vs(ii) = cdft_reg%Vs(ii) + dV

         call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
         call SCF(energ, fock_a, rho_a, fock_b, rho_b)

         
         call cdft_get_spin_constraints(cdft_reg%Cs)
         if (cdft_chrg) then
            call cdft_get_chrg_constraints(cdft_reg%Cc)
            do jj = 1, cdft_reg%n_regions
               jacob(jj,ii) = (cdft_reg%Cc(jj) - cdft_reg%Cc_old(jj)) / dV
            enddo
            do jj = 1, cdft_reg%n_regions
               jacob(jj+cdft_reg%n_regions,ii) = &
                           (cdft_reg%Cs(jj) - cdft_reg%Cs_old(jj)) / dV
            enddo
         else
            do jj = 1, cdft_reg%n_regions
               jacob(jj,ii) = (cdft_reg%Cs(jj) - cdft_reg%Cs_old(jj)) / dV
            enddo
         endif
         cdft_reg%Vs = cdft_reg%Vs_old
      enddo

      cdft_reg%Cs = cdft_reg%Cs_old
      if (cdft_chrg) cdft_reg%Cc = cdft_reg%Cc_old
   endif

   if (size(jacob,1) > 1) then
      allocate(IPIV(size(jacob,1), size(jacob,1)), WORK(1))
      jacob_i = jacob

      call DGETRI(size(jacob,1), jacob_i, size(jacob,1), IPIV, WORK, -1, INFO)
      LWORK = int(WORK(1))
      deallocate(WORK)
      allocate(WORK(LWORK))

      call DGETRF(size(jacob,1), size(jacob,1), jacob_i, size(jacob,1), IPIV, INFO)
      call DGETRI(size(jacob,1), jacob_i, size(jacob,1), IPIV, WORK, LWORK, &
                  INFO)

      deallocate(IPIV, WORK)
   endif

   print*, "JACOB" , jacob
end subroutine cdft_get_jacobian

! Propagates the constraint potentials.
! First, it gets dW/dVk (W being total energy and Vk the constrained
! potential) which equals the constraint value. Then stores the derivative
! and performs another cycle, extrapolating a new Vk from previous iterations.
subroutine cdft_set_potential()
   use cdft_data, only: cdft_chrg, cdft_spin, jacob, jacob_i, cdft_reg, j_step
   implicit none
   integer      :: ii, jj, sp_ind

   sp_ind = 0
   if (cdft_chrg .and. cdft_spin) sp_ind = cdft_reg%n_regions

   if (cdft_chrg) then
      if (size(jacob,1) > 1) then
         cdft_reg%Vc = cdft_reg%Vc_old
         do ii = 1, cdft_reg%n_regions
         do jj = 1, cdft_reg%n_regions
            cdft_reg%Vc(ii) = cdft_reg%Vc(ii) - j_step * &
                              cdft_reg%Cc_old(jj) * jacob_i(ii,jj)
         enddo
         enddo

         if (cdft_spin) then
            do ii = 1, cdft_reg%n_regions
            do jj = 1, cdft_reg%n_regions
               cdft_reg%Vc(ii) = cdft_reg%Vc(ii) - j_step * &
                                 cdft_reg%Cs(jj) * jacob_i(ii, jj+sp_ind)
            enddo
            enddo
         endif
      else
         cdft_reg%Vc(1) = cdft_reg%Vc_old(1) - j_step * &
                          cdft_reg%Cc_old(1) / jacob(1,1)
      endif
   endif

   if (cdft_spin) then
      if (size(jacob,1) > 1) then
         cdft_reg%Vs = cdft_reg%Vs_old
    
         if (cdft_chrg) then
            do ii = 1, cdft_reg%n_regions
            do jj = 1, cdft_reg%n_regions
               cdft_reg%Vs(ii) = cdft_reg%Vs(ii) - j_step * &
                                 cdft_reg%Cc(jj) * jacob_i(ii,jj)
            enddo
            enddo
         endif

         do ii = 1, cdft_reg%n_regions
            do jj = 1, cdft_reg%n_regions
               cdft_reg%Vs(ii) = cdft_reg%Vs(ii) - j_step * &
                                 cdft_reg%Cs(jj) * jacob_i(ii,jj+sp_ind)
            enddo
         enddo
      else
         cdft_reg%Vs(1) = cdft_reg%Vs_old(1) - j_step * &
                          cdft_reg%Cs_old(1) / jacob(1,1)
      endif
   endif

   call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
   print*, "CDFT Vc: ", cdft_reg%Vc
end subroutine cdft_set_potential

! Checks if CDFT converged.
subroutine cdft_check_conver(rho_new, rho_old, converged, cdft_iter)
   use cdft_data, only: cdft_reg, cdft_chrg, cdft_spin, j_step, c_prev
   implicit none
   real(kind=8), intent(in)  :: rho_new(:), rho_old(:)
   integer     , intent(in)  :: cdft_iter
   logical     , intent(out) :: converged 

   real(kind=8) :: rho_diff, const_sum
   integer      :: jj
   
   rho_diff = 0.0D0
   do jj = 1 , size(rho_new,1)
      rho_diff  = rho_diff + (rho_new(jj) - rho_old(jj)) * &
                             (rho_new(jj) - rho_old(jj))
   enddo
   rho_diff = sqrt(rho_diff) / dble(size(rho_new,1))

   const_sum = 0.0D0
   if (cdft_chrg) then
      call cdft_get_chrg_constraints(cdft_reg%Cc)
      do jj = 1, cdft_reg%n_regions
         const_sum   = const_sum + abs(cdft_reg%Cc(jj))
      enddo
      print*, "Constraints", cdft_reg%Cc
   endif
   if (cdft_spin) then
      call cdft_get_spin_constraints(cdft_reg%Cs)
      do jj = 1, cdft_reg%n_regions
         const_sum   = const_sum + abs(cdft_reg%Cs(jj))
      enddo
   endif

   write(*,*) "CDFT Iteration: ", cdft_iter, "Rho: ", rho_diff, &
              "Constraint: ", const_sum
   converged = .false.
   if ((rho_diff < 1D-4) .and. (const_sum < 1D-8)) converged = .true.

   if (const_sum > c_prev) then
      j_step = j_step * 0.8D0
   else
      j_step = j_step * 1.2D0
   endif
   c_prev = const_sum
   
end subroutine cdft_check_conver

! Adds CDFT terms to total energy.
subroutine cdft_add_energy(energ)
   use cdft_data, only: cdft_chrg, cdft_spin, cdft_reg
   implicit none
   real(kind=8), intent(inout) :: energ
   integer :: ii

   if (cdft_chrg) then
      call cdft_get_chrg_constraints(cdft_reg%Cc)
      do ii = 1, cdft_reg%n_regions
         energ   = energ + cdft_reg%Cc(ii) * cdft_reg%Vc(ii)
      enddo
   endif

   if (cdft_spin) then
      call cdft_get_spin_constraints(cdft_reg%Cs)
      do ii = 1, cdft_reg%n_regions
         energ   = energ + cdft_reg%Cs(ii) * cdft_reg%Vs(ii)
      enddo
   endif
end subroutine cdft_add_energy

! Gets Becke atomic charges and calculates its difference
! with the constrained total value.
subroutine cdft_get_chrg_constraints(constraints)
   use cdft_data, only: at_chrg, cdft_reg
   implicit none
   real(kind=8), intent(inout) :: constraints(:)
   integer                     :: c_index, region

   call g2g_get_becke_dens(at_chrg)

   constraints = 0.0D0
   do region  = 1, cdft_reg%n_regions
      do c_index = 1, cdft_reg%natom(region)
         constraints(region) = constraints(region) + at_chrg( &
                                       cdft_reg%atoms(region,c_index))
      enddo
      constraints(region) = cdft_reg%chrg(region) - constraints(region)
   enddo
end subroutine cdft_get_chrg_constraints

! Gets Becke atomic spin and calculates its difference
! with the constrained total value.
subroutine cdft_get_spin_constraints(constraints)
   use cdft_data, only: cdft_reg, at_spin
   implicit none
   real(kind=8), intent(inout) :: constraints(:)
   integer                     :: c_index, region

   call g2g_get_becke_spin(at_spin)

   constraints = 0.0D0
   do region  = 1, cdft_reg%n_regions
      do c_index = 1, cdft_reg%natom(region)
         constraints(region) = constraints(region) + at_spin( &
                                       cdft_reg%atoms(region,c_index))
      enddo
      constraints(region) = constraints(region) - cdft_reg%spin(region)
   enddo
end subroutine cdft_get_spin_constraints

end module cdft_subs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
