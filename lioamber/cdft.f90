!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% CONSTRAINED DFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains routines (cdft_subs) and variables (cdft_data) for        !
! Constrained DFT (CDFT) calculations. This implementation uses the Becke      !
! partitioning scheme for the constraints, as used by Holmberg and Laasoen     !
! (JCTC 2016, DOI: 10.1021/acs.jctc.6b01085).                                  !
!                                                                              !
! The following subroutines are called externally:                             !
!   * cdft_options_check: checks internal consistency in input options.        !
!   * cdft_input_read:    reads CDFT block input.                              !
!   * CDFT:               main CDFT loop outside SCF.                          !
!                                                                              !
! The following subroutines are only called internally:                        !
!   * cdft_get_constraints: gets sum(atomic_spin/charges) - spin/charge_target.!
!   * cdft_add_energy:      adds CDFT terms to energy.                         !
!   * cdft_initialise:      allocates arrays for CDFT.                         !
!   * cdft_finalise:        deallocates arrays for CDFT.                       !
!   * cdft_get_deltaV:      gets the increment for bias potentials.            !
!   * cdft_set_potential:   sets new potential for grid integration.           !
!   * cdft_check_conver:    checks CDFT convergence.                           !
!                                                                              !
!------------------------------------------------------------------------------!
! How to input CDFT options:                                                   !
!   CDFT options must be provided in the LIO input file (lio.in) but outside   !
!   namelists. Said file should contain the following sections:                !
!                                                                              !
!{CDFT}                                                                        !
!  N_REG CONST_CHARGE CONST_SPIN                                               !
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

   real(kind=8), allocatable :: at_chrg(:)       ! List of atomic charges.
   real(kind=8), allocatable :: at_spin(:)       ! List of atomic spin charges.
   real(kind=8), allocatable :: jacob(:,:)       ! Jacobian for advancing Vi
   integer                   :: sp_idx  = 0      ! Starting index for spin.

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
      real(kind=8), allocatable :: Vmix(:)    ! Array with both potentials.
      real(kind=8), allocatable :: cst(:)     ! Charge/Spin constraints.
      
      ! Arrays for Jacobian calculation and propagation
      real(kind=8), allocatable :: Vc_old(:)  ! Charge potentials.
      real(kind=8), allocatable :: Vs_old(:)  ! Spin potentials.
      real(kind=8), allocatable :: Vm_old(:)  ! Array with both potentials.
      real(kind=8), allocatable :: cst_old(:) ! Charge/Spin constraint value.

   end type cdft_region_data

   type(cdft_region_data) :: cdft_reg
end module cdft_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module cdft_subs
   implicit none
   private
   public :: cdft_input_read
   public :: cdft_options_check
   public :: CDFT
contains

!%%%% PUBLIC SUBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
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

   write(*,'(A)') "CDFT input found, reading options."
   ! Starts reading CDFT data.
   read(input_UID,*) cdft_reg%n_regions, inp_chrg, inp_spin
   if (inp_chrg == 1) cdft_chrg = .true.
   if (inp_spin == 1) cdft_spin = .true.
   write(*,'(A21,I3,A21,L2,A19,L2)')"  Number of regions: ",cdft_reg%n_regions,&
                                    " | Constrain charge: ", cdft_chrg, &
                                    " | Constrain spin: ", cdft_spin

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

! Checks if doing CDFT and sets energy_all_iterations to true in order to
! also calculate Becke charges in each iteration step.
subroutine cdft_options_check(do_becke, open_shell)
   use cdft_data, only: cdft_chrg, cdft_spin, doing_cdft

   implicit none
   logical, intent(inout) :: do_becke, open_shell

   if ((cdft_chrg) .or. (cdft_spin)) then
      doing_cdft     = .true.
      do_becke       = .true.
   endif

   if (cdft_spin) open_shell = .true.
end subroutine cdft_options_check

! Performs the CDFT iterative procedure.
subroutine CDFT(fock_a, rho_a, fock_b, rho_b, Pmat_vec, natom)
   use typedef_operator, only: operator
   use converger_data  , only: told

   implicit none
   integer       , intent(in)    :: natom
   real(kind=8)  , intent(inout) :: Pmat_vec(:)
   type(operator), intent(inout) :: fock_a, rho_a, fock_b, rho_b

   integer      :: cdft_iter, max_cdft_iter
   logical      :: cdft_converged = .false.
   real(kind=8) :: energ
   real(kind=8), allocatable :: Pmat_old(:)

   max_cdft_iter = 100
   cdft_iter     = 0
   allocate(Pmat_old(size(Pmat_vec,1)))
   call cdft_initialise(natom)
   do while ((.not. cdft_converged) .and. (cdft_iter < max_cdft_iter))
      cdft_iter = cdft_iter +1
      Pmat_old  = Pmat_vec
      call SCF(energ, fock_a, rho_a, fock_b, rho_b)
      call cdft_check_conver(Pmat_vec, Pmat_old, cdft_converged, &
                             cdft_iter, energ, told)

      if (.not. cdft_converged) then
         ! Calculates perturbations and Jacobian.
         call cdft_get_deltaV(fock_a, rho_a, fock_b, rho_b)
         call cdft_set_potential()
      endif
   enddo

   call cdft_finalise()
   deallocate(Pmat_old)
end subroutine CDFT

!%%%% PRIVATE SUBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Initializes arrays.
subroutine cdft_initialise(n_atoms)
   use cdft_data, only: cdft_reg, cdft_chrg, cdft_spin, &
                        at_chrg, at_spin, jacob, sp_idx
   implicit none
   integer, intent(in) :: n_atoms
   integer             :: J_size = 0

   if (cdft_chrg .and. cdft_spin) then
      J_size = 2 * cdft_reg%n_regions
      sp_idx = cdft_reg%n_regions
   else if (cdft_chrg .or. cdft_spin) then
      J_size = cdft_reg%n_regions
      sp_idx = 0
   endif

   if (cdft_chrg .or. cdft_spin) then
      if (allocated(at_chrg))         deallocate(at_chrg)
      if (allocated(cdft_reg%Vc))     deallocate(cdft_reg%Vc)
      if (allocated(cdft_reg%Vc_old)) deallocate(cdft_reg%Vc_old)
      allocate(at_chrg(n_atoms))
      allocate(cdft_reg%Vc(cdft_reg%n_regions))
      allocate(cdft_reg%Vc_old(cdft_reg%n_regions))

      if (allocated(at_spin))         deallocate(at_spin)
      if (allocated(cdft_reg%Vs))     deallocate(cdft_reg%Vs)
      if (allocated(cdft_reg%Vs_old)) deallocate(cdft_reg%Vs_old)
      allocate(at_spin(n_atoms))  
      allocate(cdft_reg%Vs(cdft_reg%n_regions))
      allocate(cdft_reg%Vs_old(cdft_reg%n_regions))

      if (allocated(cdft_reg%cst))     deallocate(cdft_reg%cst)
      if (allocated(cdft_reg%cst_old)) deallocate(cdft_reg%cst_old)
      if (allocated(cdft_reg%Vmix))    deallocate(cdft_reg%Vmix)
      if (allocated(cdft_reg%Vm_old))  deallocate(cdft_reg%Vm_old)
      allocate(cdft_reg%cst(J_size))
      allocate(cdft_reg%cst_old(J_size))
      allocate(cdft_reg%Vmix(J_size))
      allocate(cdft_reg%Vm_old(J_size))
   endif

   if (allocated(jacob)) deallocate(jacob)
   allocate(jacob(J_size, J_size))

   call g2g_cdft_init(cdft_chrg, cdft_spin, cdft_reg%n_regions, &
                      cdft_reg%max_nat, cdft_reg%natom, cdft_reg%atoms)
   cdft_reg%Vc = 0.0D0                      
   cdft_reg%Vs = 0.0D0
   call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
end subroutine cdft_initialise

! Deallocates arrays.
subroutine cdft_finalise()
   use cdft_data, only: cdft_reg, at_spin, at_chrg, jacob
   implicit none

   if (allocated(jacob))            deallocate(jacob)
   if (allocated(at_chrg))          deallocate(at_chrg)
   if (allocated(at_spin))          deallocate(at_spin)
   if (allocated(cdft_reg%Vc))      deallocate(cdft_reg%Vc)
   if (allocated(cdft_reg%Vs))      deallocate(cdft_reg%Vs)
   if (allocated(cdft_reg%natom))   deallocate(cdft_reg%natom)
   if (allocated(cdft_reg%atoms))   deallocate(cdft_reg%atoms)
   if (allocated(cdft_reg%chrg))    deallocate(cdft_reg%chrg)
   if (allocated(cdft_reg%spin))    deallocate(cdft_reg%spin)
   if (allocated(cdft_reg%Vc_old))  deallocate(cdft_reg%Vc_old)
   if (allocated(cdft_reg%Vs_old))  deallocate(cdft_reg%Vs_old)
   if (allocated(cdft_reg%Vmix))    deallocate(cdft_reg%Vmix)
   if (allocated(cdft_reg%Vm_old))  deallocate(cdft_reg%Vm_old)
   if (allocated(cdft_reg%cst))     deallocate(cdft_reg%cst)
   if (allocated(cdft_reg%cst_old)) deallocate(cdft_reg%cst_old)

   call g2g_cdft_finalise()
end subroutine cdft_finalise

! Gets the jacobian matrix by making a small perturbation in each direction.
! This is done in order to propagate the constraint potentials Vc and Vs
! by means of Newton's method. Instead of calculating J-1, we solve an
! alternative problem resulting in ΔVi (i=c,s). 
subroutine cdft_get_deltaV(fock_a, rho_a, fock_b, rho_b)
   use typedef_operator, only: operator
   use cdft_data       , only: cdft_chrg, cdft_spin, cdft_reg, &
                               jacob, sp_idx
   
   implicit none
   type(operator), intent(inout) :: fock_a, rho_a, fock_b, rho_b
   integer      :: ii, jj
   real(kind=8) :: energ, dV

   ! Variables for LAPACK
   integer                   :: LWORK, INFO
   real(kind=8), allocatable :: WORK(:)

   jacob = 0.0D0
   call cdft_get_constraints()
   cdft_reg%cst_old = cdft_reg%cst

   if (cdft_chrg) then
      cdft_reg%Vc_old = cdft_reg%Vc

      do ii = 1, cdft_reg%n_regions
         dV = cdft_reg%cst(ii)
         if (dV > 0.01D0) dV = 0.01D0
         cdft_reg%Vc(ii) = cdft_reg%Vc(ii) + dV

         call g2g_cdft_set_v(cdft_reg%Vc, cdft_reg%Vs)
         call SCF(energ, fock_a, rho_a, fock_b, rho_b)

         call cdft_get_constraints()
         do jj = 1, cdft_reg%n_regions
            jacob(jj,ii) = (cdft_reg%cst(jj) - cdft_reg%cst_old(jj)) / dV
         enddo
         if (cdft_spin) then
            do jj = 1+sp_idx, cdft_reg%n_regions+sp_idx
               jacob(jj,ii) = (cdft_reg%cst(jj) - cdft_reg%cst_old(jj)) / dV
            enddo
         endif
         cdft_reg%Vc = cdft_reg%Vc_old
      enddo
   endif

   if (cdft_spin) then
      cdft_reg%Vs_old = cdft_reg%Vs

      do ii = 1, cdft_reg%n_regions
         dV = cdft_reg%cst(ii+sp_idx)
         if (dV > 0.01D0) dV = 0.01D0
         cdft_reg%Vs(ii) = cdft_reg%Vs(ii) + dV

         call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
         call SCF(energ, fock_a, rho_a, fock_b, rho_b)
         
         call cdft_get_constraints()
         if (cdft_chrg) then
            do jj = 1, cdft_reg%n_regions
               jacob(jj,ii+sp_idx) = (cdft_reg%cst(jj) - &
                                      cdft_reg%cst_old(jj)) / dV
            enddo
         endif
         
         do jj = 1+sp_idx, cdft_reg%n_regions+sp_idx
            jacob(jj,ii+sp_idx) = (cdft_reg%cst(jj) - &
                                   cdft_reg%cst_old(jj)) / dV
         enddo
         cdft_reg%Vs = cdft_reg%Vs_old
      enddo
   endif
   cdft_reg%cst = cdft_reg%cst_old

   ! Alternative to invert J-1: Solving A*x = B where A is the jacobian
   ! and B is the negative constraints array. The result is an array containing
   ! Xn+1 - Xn for each constraint. This is totally equivalent to solve
   ! Xn+1 = Xn - J^(-1) * Cst, with A = J and B = -Cst, but much less costly.
   if (size(jacob,1) > 1) then
      cdft_reg%Vmix = -cdft_reg%cst 
      allocate(WORK(1))
      call dgels('N', size(jacob,1), size(jacob,1), 1, jacob, size(jacob,1), &
                 cdft_reg%Vmix, size(jacob,1), WORK, -1, INFO)
      LWORK = int(WORK(1))
      deallocate(WORK)
      allocate(WORK(LWORK))
      call dgels('N', size(jacob,1), size(jacob,1), 1, jacob, size(jacob,1), &
                 cdft_reg%Vmix, size(jacob,1), WORK, LWORK, INFO)
      deallocate(WORK)
   endif
end subroutine cdft_get_deltaV

! Propagates the constraint potentials by means of Newton's method. Vmix,
! which contains ΔVc and ΔVs, is obtained in the previous routine.
subroutine cdft_set_potential()
   use cdft_data, only: cdft_chrg, cdft_spin, jacob, cdft_reg, sp_idx
   implicit none

   if (size(jacob,1) == 1) then
      if (cdft_chrg) then
         cdft_reg%Vc(1) = cdft_reg%Vc_old(1) - cdft_reg%cst(1) / jacob(1,1)
      else if (cdft_spin) then
         cdft_reg%Vs(1) = cdft_reg%Vs_old(1) - cdft_reg%cst(1) / jacob(1,1)
      endif
   else
      if (cdft_chrg) cdft_reg%Vm_old(1:size(cdft_reg%Vc,1)) = &
                     cdft_reg%Vc_old(:)
      if (cdft_spin) &
               cdft_reg%Vm_old((sp_idx+1):(sp_idx+size(cdft_reg%Vs,1))) = &
               cdft_reg%Vs_old(:)

      cdft_reg%Vmix = cdft_reg%Vmix + cdft_reg%Vm_old

      if (cdft_chrg) cdft_reg%Vc = cdft_reg%Vmix(1:size(cdft_reg%Vc,1))
      if (cdft_spin) cdft_reg%Vs = &
                     cdft_reg%Vmix((sp_idx+1):(sp_idx+size(cdft_reg%Vs,1)))
   endif

   call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
end subroutine cdft_set_potential

! Checks if CDFT converged.
subroutine cdft_check_conver(rho_new, rho_old, converged, cdft_iter, ener, &
                             rho_crit)
   use cdft_data, only: cdft_reg
   implicit none
   real(kind=8), intent(in)    :: rho_new(:), rho_old(:), rho_crit
   integer     , intent(in)    :: cdft_iter
   logical     , intent(out)   :: converged
   real(kind=8), intent(inout) :: ener

   real(kind=8) :: rho_diff, c_max
   integer      :: jj
   
   rho_diff = 0.0D0
   do jj = 1 , size(rho_new,1)
      rho_diff  = rho_diff + (rho_new(jj) - rho_old(jj)) * &
                             (rho_new(jj) - rho_old(jj))
   enddo
   rho_diff = sqrt(rho_diff) / dble(size(rho_new,1))

   call cdft_get_constraints()
   c_max = maxval(abs(cdft_reg%cst))

   call cdft_add_energy(ener)
   write(*,'(A)') "CDFT Convergence status:" 
   write(*,*) "Iteration n°:      ", cdft_iter
   write(*,*) "Energy:            ", ener
   write(*,*) "ΔRho:              ", rho_diff
   write(*,*) "Constraint values: ", cdft_reg%cst
   write(*,*) "Charge potential:  ", cdft_reg%Vc
   write(*,*) "Spin potential:    ", cdft_reg%Vs
   converged = .false.
   if ((rho_diff < rho_crit) .and. (c_max < 1D-5)) converged = .true.
end subroutine cdft_check_conver

! Adds CDFT terms to total energy.
subroutine cdft_add_energy(energ)
   use cdft_data, only: cdft_chrg, cdft_spin, cdft_reg, sp_idx
   implicit none
   real(kind=8), intent(inout) :: energ
   integer :: ii

   call cdft_get_constraints()
   if (cdft_chrg) then
      do ii = 1, cdft_reg%n_regions
         energ = energ + cdft_reg%cst(ii) * cdft_reg%Vc(ii)
      enddo
   endif

   if (cdft_spin) then
      do ii = 1, cdft_reg%n_regions
         energ = energ + cdft_reg%cst(ii+sp_idx) * cdft_reg%Vs(ii)
      enddo
   endif
end subroutine cdft_add_energy

! Gets Becke atomic spin and calculates its difference
! with the constrained total value.
subroutine cdft_get_constraints()
   use cdft_data, only: cdft_reg, at_spin, at_chrg, &
                        cdft_spin, cdft_chrg, sp_idx
   implicit none
   integer :: c_index, region

   cdft_reg%cst = 0.0D0
   if (cdft_chrg) then
      call g2g_get_becke_dens(at_chrg)

      do region = 1, cdft_reg%n_regions
         do c_index = 1, cdft_reg%natom(region)
            cdft_reg%cst(region) = cdft_reg%cst(region) + &
                                   at_chrg(cdft_reg%atoms(region,c_index))
         enddo
         cdft_reg%cst(region) = cdft_reg%chrg(region) - cdft_reg%cst(region)
      enddo
   endif

   if (cdft_spin) then
      call g2g_get_becke_spin(at_spin)

      do region = 1, cdft_reg%n_regions
         do c_index = 1, cdft_reg%natom(region)
            cdft_reg%cst(region+sp_idx)= cdft_reg%cst(region+sp_idx) + &
                                         at_spin(cdft_reg%atoms(region,c_index))
         enddo
         cdft_reg%cst(region+sp_idx) = cdft_reg%cst(region+sp_idx) - &
                                       cdft_reg%spin(region)
                                       
      enddo
   endif
end subroutine cdft_get_constraints

end module cdft_subs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
