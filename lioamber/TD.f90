!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% REAL TIME-TDDFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! 1992 - Dario Estrin (source code)                                            !
! 2012 - Nano, Dario, Uriel, Damian (main routine, external fields)            !
! 2015 - Fran (magnus)                                                         !
! 2017 - Charly, Gonzalo, Fede (reformating)                                   !
!                                                                              !
! This subroutine takes the converged density matrix from an SCF calculation   !
! and evolves it in time. In the input file the total number of propagation    !
! steps is specified (nstep) as well as the time of each evolution step        !
! (tdstep).                                                                    !
! This implementation has two alternatives to evolve the density in time. The  !
! first one (propagator=1) is the Verlet algorithm that uses a convination of  !
! Liouville von Newmann expresion for the time derivative of the density matrix!
! and a first order Taylor expansion of the density matrix. The second one     !
! (propagator=2) is the Magnus propagation scheme that uses Backer Campbell    !
! Hausdorff (BCH) formula. For this reason when Magnus is used the number of   !
! total conmutators in the BCH espansion has to be specified (NBCH,            !
! default=10).                                                                 !
! A narrow gaussian type electric field can be introduced during the time      !
! evolution in order to excite all electronic frequencies with the same        !
! intensity. Once this perturbation is turned on (Field=t, exter=t) each       !
! component of the external electric field has to be specified in the input    !
! file (Fx,Fy,Fz).                                                             !
! In each step of the propagation the cartesian components of the sistems      !
! dipole are stored in files x.dip, y.dip, z.dip.                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module td_data
   implicit none
   integer :: td_rst_freq = 500
   integer :: timedep   = 0
   integer :: ntdstep   = 0
   real*8  :: tdstep    = 2.0D-3
   logical :: tdrestart = .false.
   logical :: writedens = .false.
   real*8  :: pert_time = 2.0D-1
end module td_data

module time_dependent
   implicit none
contains

subroutine TD()
   use garcha_mod    , only: M, Md, NBCH, propagator, RMM, NCO, Iz, igrid2, r, &
                             Nuc, nsol, pc, X, Smat, MEMO, sol, natom, sqsm,   &
                             Nunp, ntatom, ncont, nshell
   use td_data       , only: td_rst_freq, tdstep, ntdstep, tdrestart, &
                             writedens, pert_time
   use field_data    , only: field, fx, fy, fz
   use field_subs    , only: field_setup_old, field_finalize
   use transport_data, only: transport_calc
   use transport_subs, only: transport_rho_trace, transport_generate_rho,      &
                             transport_init, transport_population
   use dftb_data     , only: dftb_calc, MDFTB, MTB
   use dftb_subs     , only: dftb_td_init, dftb_output
   use fileio        , only: write_td_restart_verlet, write_td_restart_magnus, &
                             read_td_restart_verlet , read_td_restart_magnus
#ifdef CUBLAS
   use cublasmath    , only: basechange_cublas
#endif

   implicit none
   real*8  :: E, En, E1, E2, E1s, Es, Ens = 0.0D0, Ex, t, dt_magnus, dt_lpfrg
   integer :: MM, MMd, M2, M5, M13, M15, M11, LWORK, igpu, info, istep, icount,&
              jcount
   integer :: lpfrg_steps = 200, chkpntF1a = 185, chkpntF1b = 195
   logical :: is_lpfrg
   character(len=20) :: restart_filename

   real*8 , allocatable, dimension(:)   :: factorial, WORK
   real*8 , allocatable, dimension(:,:) :: Xmm, Xmat,Xtrans, Ytrans, fock, F1a,&
                                           F1b, overlap, Ymat
   real*8 , allocatable, dimension(:,:) :: fock_0

#ifdef TD_SIMPLE
   complex*8, allocatable, dimension(:,:) :: rho_0
#else
   complex*16, allocatable, dimension(:,:) :: rho_0
#endif
! Precision options.
#ifdef TD_SIMPLE
   complex*8  :: Im = (0.0E0,2.0E0)
   complex*8 , allocatable, dimension(:,:) :: rho, rho_aux, rhonew, rhold
#else
   complex*16 :: Im = (0.0D0,2.0D0)
   complex*16, allocatable, dimension(:,:) :: rho, rho_aux, rhonew, rhold
#endif

! CUBLAS options.
#ifdef CUBLAS
   integer   :: sizeof_real, sizeof_complex
   integer*8 :: devPtrX, devPtrY, devPtrXc
   parameter(sizeof_real    = 8)
#ifdef TD_SIMPLE
   parameter(sizeof_complex = 8)
#else
   parameter(sizeof_complex = 16)
#endif
#endif
!DFTB: M_in controls de size of the bigest matrices for DFTB, ii and jj are only
!counters, and traza is for the control of the trace of density matrix
   integer :: M_in, ii,jj

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% TD INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   ! RMM initializations
   ! M3 - Fold, M5 - S and then F, M7 - G, M9 - Gm, M11 - H
   ! M13 - W (matrix eigenvalues),  M15 - Aux vector for ESSL. M17 - Used
   ! in least squares, M18 - MO vectors (deprecated), M19 - weights (optional),
   ! M20 - 2e integrals if MEMO=t
   MM  = M *(M+1) /2     ; MMd = Md*(Md+1)/2
   M2  = 2*M             ; M5  = 1 + 2*MM
   M11 = M5 + MM + 2*MMd ; M13 = M11 + MM
   M15 = M13 + M

   call g2g_timer_start('TD')
   call g2g_timer_start('td-inicio')
   open(unit = 134, file = "dipole_moment_td")

!------------------------------------------------------------------------------!
!DFTB: defining DFTB matrix size

      if (dftb_calc) then
         M_in=MDFTB
      else
         M_in=M
      end if
!------------------------------------------------------------------------------!

   ! Checks and performs allocations.
   call td_allocate_all(M_in, M, NBCH, propagator, F1a, F1b, fock, rho, rho_aux, &
                        rhold, rhonew, rho_0, fock_0, sqsm, Xmm, Xmat, Xtrans,   &
                        Ymat, Ytrans, factorial)
#ifdef CUBLAS
   call td_allocate_cublas(M_in, sizeof_real, devPtrX, devPtrY)
#endif

   ! Initialises propagator-related parameters and other variables.
   call td_initialise(propagator, tdstep, NBCH, dt_lpfrg, dt_magnus, factorial,&
                      NCO, Nunp, natom, Iz)

   ! TD restart reading.
   if (tdrestart) then
      restart_filename = 'td_in.restart'
      if (propagator.eq.2) then
         call read_td_restart_magnus(rho, F1a, F1b, M_in, restart_filename)
      else
         call read_td_restart_verlet(rho, M_in, restart_filename)
      endif
      call sprepack_ctr('L', M, RMM, rho_0)
   else
      ! Read the density matrix stored in RMM(1,2,3,...,MM) into rho matrix.
      call spunpack_rtc('L', M, RMM, rho_0)
   endif

!------------------------------------------------------------------------------!
! DFTB: Initialize DFTB variables for TD
   if (dftb_calc) then
      call dftb_td_init (M , rho, rho_0, overlap, RMM(M5))
   else
      rho=rho_0
   end if
!------------------------------------------------------------------------------!

   ! Proper TD calculation start.
   write(*,*) 'Starting TD calculation...'
   ! Create integration grid for XC, assigning points to groups (spheres/cubes)
   ! and significant functions to groups, also calculating point weights.
   if (field) call field_setup_old(pert_time, 1, fx, fy, fz)
   call td_integration_setup(igrid2, igpu)
   call td_integral_1e(E1, En, E1s, Ens, MM, igpu, nsol, RMM, RMM(M11), r, pc, &
                       ntatom)
   ! Initialises transport if required.
   if (transport_calc) call transport_init(M, natom, Nuc, RMM(M5), overlap, rho)

   ! Diagonalizes Smat, stores transformation matrix Xmm and calculates the
   ! transposed matrices Xtrans and Ytrans
   call td_overlap_diag(M_in, M, Smat, RMM(M13), X, Xmat ,Xtrans, Ymat, Ytrans,&
                        Xmm)
   ! Here rho_aux is used as an auxiliar matrix for CUBLAS. Then we need to
   ! restore its original value. Rho is then transformed to the orthonormal
   ! basis with either matmul intrinsic or CUBLAS.
#ifdef CUBLAS
   call td_bc_rho_cu(M_in, sizeof_real, sizeof_complex, devPtrX, devPtrY, &
                     devPtrXc, Xmat, Ymat, rho, rho_aux)
#else
   call td_bc_rho(M_in, rho, Ymat, Ytrans)
#endif
   ! Precalculate three-index (two in MO basis, one in density basis) matrix
   ! used in density fitting /Coulomb F element calculation here (t_i in Dunlap)
   call td_coulomb_precalc(igpu, MEMO)

   call g2g_timer_stop('td-inicio')
   ! End of TD initialization.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% TD EVOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   do 999 istep = 1, ntdstep
      call g2g_timer_start('TD step')
      call g2g_timer_sum_start("TD - TD Step")
      ! Checks if step is a leapfrog step.
      call td_check_prop(is_lpfrg, propagator, istep, lpfrg_steps, tdrestart)
      call td_get_time(t, tdstep, istep, propagator, is_lpfrg)

      call g2g_timer_sum_start("TD - TD Step Energy")
      call td_calc_energy(E, E1, E2, En, Ex, Es, MM, RMM, RMM(M11), is_lpfrg,  &
                          transport_calc, sol, t/0.024190D0)
      call g2g_timer_sum_pause("TD - TD Step Energy")
      write(*,*) 'TD - Step: ', istep, " Energy : ", E

      ! Verlet or Magnus Propagation
      ! In Verlet, we obtain the Fock matrix in the molecular orbital (MO)
      ! basis, where U matrix with eigenvectors of S, and s is vector with
      ! eigenvalues. In the first step of the propagation we extrapolate rho
      ! back in time using Verlet algorithm to calculate rhold.
      ! Both propagations include transport.

      ! After propagation, we transform the density to the atomic orbital basis
      ! and take the real part of it. The imaginary part of the density can be
      ! discarded since for a basis set of purely real functions the fock matrix
      ! is real and symetric, and depends only on the real part of the complex
      ! density matrix. (This will not be true in the case of hybrid
      ! functionals)

      call g2g_timer_sum_start("TD - Propagation")
#ifdef CUBLAS
      if (is_lpfrg) then
         call td_bc_fock_cu(M_in,M, MM, RMM(M5), fock, devPtrX,natom, nshell,ncont, istep)

         call td_verlet_cu(M, M_in, fock, rhold, rho, rhonew, istep, Im, dt_lpfrg, &
                           transport_calc, natom, Nuc, Iz, overlap, sqsm,    &
                           rho_aux, devPtrY,devPtrXc)


         if (propagator.eq.2) then
            if (istep.eq.chkpntF1a) F1a = fock
            if (istep.eq.chkpntF1b) F1b = fock
         endif
      else


         call td_magnus_cu(M, fock, F1a, F1b, rho, rhonew, devPtrX, devPtrXc,  &
                           factorial, NBCH, dt_magnus, natom, transport_calc,  &
                           Nuc, Iz, istep, overlap, sqsm, rho_aux, devPtrY,    &
                           t/0.024190D0, M_in, nshell, ncont)
      endif

      call g2g_timer_start('complex_rho_on_to_ao-cu')
      rho_aux = basechange_cublas(M_in, rho, devPtrXc, 'inv')
      call g2g_timer_stop('complex_rho_on_to_ao-cu')

#else
      if (is_lpfrg) then
         call td_bc_fock(M_in,M, MM, RMM(M5), fock, Xtrans, Xmm, natom, nshell,ncont, istep)
         call td_verlet(M, M_in, fock, rhold, rho, rhonew, istep, Im, dt_lpfrg,   &
                        transport_calc, natom, Nuc, Iz, overlap,sqsm, rho_aux,    &
                        Xmat, Xtrans)
         if (propagator.eq.2) then
            if (istep.eq.chkpntF1a) F1a = fock
            if (istep.eq.chkpntF1b) F1b = fock
         endif
      else
         call td_magnus(M, fock, F1a, F1b, rho, rhonew, factorial, NBCH,  &
                        dt_magnus, natom, transport_calc, Nuc, Iz, istep, &
                        overlap, sqsm, rho_aux, X, Xtrans, t/0.024190D0, M_in,&
                        nshell, ncont)
      endif

      call g2g_timer_start('complex_rho_on_to_ao')
      rho_aux = matmul(Xmat, rho)
      rho_aux = matmul(rho_aux, Xtrans)
      call g2g_timer_stop('complex_rho_on_to_ao')
#endif
      call g2g_timer_sum_pause("TD - Propagation")

      ! The real part of the density matrix in the atomic orbital basis is
      ! copied in RMM(1,2,3,...,MM) to compute the corresponding Fock matrix.
      call sprepack_ctr('L', M, RMM, rho_aux(MTB+1:MTB+M,MTB+1:MTB+M))

!DFTB: temporary DFTB don't store restarts
      ! Stores the density matrix each 500 steps as a restart.
      if ((writedens) .and. ( (mod(istep, td_rst_freq) == 0) .or. &
         (istep.eq.ntdstep) )) then
         restart_filename='td.restart'
         if (istep.eq.ntdstep) restart_filename='td_last.restart'
         if (propagator.eq.2) then
            call write_td_restart_magnus(rho, F1a, F1b, M, restart_filename)
         else
            call write_td_restart_verlet(rho, M, restart_filename)
         endif
      endif

      ! Compute the trace of the density matrix for population analysis.
      if (transport_calc) call transport_rho_trace(M, rho)

      ! Dipole Moment calculation.
      call td_dipole(t, tdstep, Fx, Fy, Fz, istep, propagator, is_lpfrg, 134)

      ! Population analysis.
      if (transport_calc) call transport_population(M, natom, Nuc, Iz, rho_aux,&
                                                    overlap, sqsm, propagator, &
                                                    is_lpfrg, istep)
      ! TD step finalization.

      if (dftb_calc) call dftb_output(M, rho_aux, overlap, istep, Iz, natom, Nuc)

      call g2g_timer_stop('TD step')
      call g2g_timer_sum_pause("TD - TD Step")

 999  continue

   ! Finalization.
   call WriteEnergies(E1, E2, En, Ens, 0, Ex, .false., 0)

#ifdef CUBLAS
   call td_finalise_cublas(devPtrX, devPtrY, devPtrXc)
#endif
   call field_finalize()

   close(134)
   call g2g_timer_stop('TD')

   return
end subroutine TD
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Subroutines used in initialisation.

subroutine td_allocate_all(M_in,M, NBCH, propagator, F1a, F1b, fock, rho, rho_aux, &
                           rhold, rhonew, rho_0, fock_0,sqsm, Xmm, Xmat, Xtrans,   &
                           Ymat, Ytrans, factorial)
   implicit none
   integer, intent(in) :: M, NBCH, propagator, M_in
   real*8, allocatable, intent(inout) :: F1a(:,:), F1b(:,:), fock(:,:),        &
                                         sqsm(:,:), Xmm(:,:), Xmat(:,:),      &
                                         Xtrans(:,:), Ymat(:,:), Ytrans(:,:),  &
                                         factorial(:), fock_0(:,:)
#ifdef TD_SIMPLE
   complex*8 , allocatable, intent(inout) :: rho(:,:), rho_aux(:,:),    &
                                             rhold(:,:), rhonew(:,:),rho_0(:,:)
#else
   complex*16, allocatable, intent(inout) :: rho(:,:), rho_aux(:,:),    &
                                             rhold(:,:), rhonew(:,:),rho_0(:,:)
#endif

   if ( allocated(fock)     ) deallocate(fock)
   if ( allocated(rho)      ) deallocate(rho)
   if ( allocated(rho_aux)  ) deallocate(rho_aux)
   if ( allocated(rhold)    ) deallocate(rhold)
   if ( allocated(rhonew)   ) deallocate(rhonew)
   if ( allocated(sqsm)     ) deallocate(sqsm)
   if ( allocated(Xmm)      ) deallocate(Xmm)
   if ( allocated(Xmat)     ) deallocate(Xmat)
   if ( allocated(Xtrans)   ) deallocate(Xtrans)
   if ( allocated(Ytrans)   ) deallocate(Ytrans)
   if ( allocated(Ymat)     ) deallocate(Ymat)
   if ( allocated(factorial)) deallocate(factorial)
   if ( allocated(fock_0)   ) deallocate(fock_0)
   if ( allocated(rho_0)    ) deallocate(rho_0)

   allocate(fock(M_in,M_in)  , rho(M_in,M_in)   , rho_aux(M_in,M_in),          &
            rhold(M_in,M_in), rhonew(M_in,M_in), sqsm(M,M), Xmm(M_in,M_in),    &
            Xtrans(M_in,M_in), Xmat(M_in,M_in), Ymat(M_in,M_in),               &
            Ytrans(M_in,M_in), factorial(NBCH), fock_0(M,M),rho_0(M,M))
   if (propagator.eq.2) then
      if ( allocated(F1a) ) deallocate(F1a)
      if ( allocated(F1b) ) deallocate(F1b)
      allocate (F1a(M_in,M_in), F1b(M_in,M_in))
   endif

   return
end subroutine td_allocate_all

subroutine td_deallocate_all(F1a, F1b, fock, rho, rho_aux, rhold, &
                             rhonew, Xmm, Xtrans, Ymat, Ytrans,factorial)
   implicit none
   real*8, allocatable, intent(inout) :: F1a(:,:), F1b(:,:), fock(:,:),    &
                                         Xmm(:,:), Xtrans(:,:), Ymat(:,:), &
                                         Ytrans(:,:), factorial(:)
#ifdef TD_SIMPLE
   complex*8 , allocatable, intent(inout) :: rho(:,:), rho_aux(:,:),    &
                                             rhold(:,:), rhonew(:,:)
#else
   complex*16, allocatable, intent(inout) :: rho(:,:), rho_aux(:,:),    &
                                             rhold(:,:), rhonew(:,:)
#endif
   if ( allocated(fock)     ) deallocate(fock)
   if ( allocated(rho)      ) deallocate(rho)
   if ( allocated(rho_aux)  ) deallocate(rho_aux)
   if ( allocated(rhold)    ) deallocate(rhold)
   if ( allocated(rhonew)   ) deallocate(rhonew)
   if ( allocated(Xmm)      ) deallocate(Xmm)
   if ( allocated(Xtrans)   ) deallocate(Xtrans)
   if ( allocated(Ytrans)   ) deallocate(Ytrans)
   if ( allocated(Ymat)     ) deallocate(Ymat)
   if ( allocated(factorial)) deallocate(factorial)
   if ( allocated(F1a)      ) deallocate(F1a)
   if ( allocated(F1b)      ) deallocate(F1b)

   return
end subroutine td_deallocate_all

subroutine td_initialise(propagator, tdstep, NBCH, dt_lpfrg, dt_magnus,        &
                         factorial, NCO, Nunp, natom, Iz)
   implicit none
   integer, intent(in)  :: propagator, NBCH, NCO, Nunp, natom, Iz(natom)
   real*8 , intent(in)  :: tdstep
   real*8 , intent(out) :: dt_lpfrg, dt_magnus, factorial(NBCH)
   integer :: icount, Nel
   real*8  :: Qc

   select case (propagator)
   ! Initialises propagator-related parameters.
      case (1)
         dt_lpfrg = tdstep
      case (2)
         dt_magnus    = tdstep
         dt_lpfrg     = tdstep*0.10D0
         factorial(1) = 1.0D0

         do icount = 2, NBCH
#ifdef CUBLAS
            factorial(icount) = 1.0D0 / icount
#else
            factorial(icount) = factorial(icount - 1) / icount
#endif
         enddo
      case default
         write(*,*) "ERROR - TD: Wrong value for propagator (td_initialise)."
         stop
   end select

   ! Calculate total atomic charge.
   Nel = 2*NCO + Nunp
   Qc  = 0.0D0

   do icount = 1, natom
      Qc = Qc + Iz(icount)
   enddo

   return
end subroutine td_initialise

subroutine td_integration_setup(igrid2, igpu)
   implicit none
   integer, intent(in)  :: igrid2
   integer, intent(out) :: igpu

   call g2g_timer_sum_start('TD - Exchange-correlation grid setup')
   call g2g_reload_atom_positions(igrid2)
   call g2g_timer_sum_stop('TD - Exchange-correlation grid setup')

   call aint_query_gpu_level(igpu)
   if (igpu.gt.1) call aint_new_step()

   return
end subroutine td_integration_setup

subroutine td_integral_1e(E1, En, E1s, Ens, MM, igpu, nsol, RMM, RMM11, r, pc, &
                          ntatom)
   use faint_cpu77, only: int1, intsol
   use mask_ecp   , only: ECP_fock
   implicit none
   integer, intent(in)    :: MM, igpu, nsol, ntatom
   real*8 , intent(in)    :: r(ntatom), pc(ntatom)
   real*8 , intent(inout) :: RMM(MM), RMM11(MM), E1, En, E1s, Ens
   integer :: icount

   E1 = 0.0D0 ; En = 0.0D0
   call g2g_timer_sum_start('TD - 1-e Fock')
   call g2g_timer_sum_start('TD - Nuclear attraction')
   call int1(En)

   call ECP_fock(MM, RMM11)
   call g2g_timer_sum_stop('TD - Nuclear attraction')

   ! 1e terms - QMMM terms.
   if ((nsol.gt.0) .or. (igpu.ge.4)) then
      call g2g_timer_sum_start('TD - QM/MM')
      if (igpu.le.1) then
         call g2g_timer_start('intsol')
         call intsol(E1s, Ens, .true.)
         call g2g_timer_stop('intsol')
      else
         call aint_qmmm_init(nsol, r, pc)
         call g2g_timer_start('aint_qmmm_fock')
         call aint_qmmm_fock(E1s, Ens)
         call g2g_timer_stop('aint_qmmm_fock')
      endif
         call g2g_timer_sum_stop('TD - QM/MM')
   endif

   E1=0.D0
   do icount = 1, MM
      E1 = E1 + RMM(icount) * RMM11(icount)
   enddo
   call g2g_timer_sum_stop('TD - 1-e Fock')
   return
end subroutine td_integral_1e

subroutine td_overlap_diag(M_in, M, Smat, eigenvalues, X_min, Xmat, Xtrans, Ymat, &
                           Ytrans, Xmm)

   use dftb_data, only:dftb_calc, MTB
   use dftb_subs, only:getXY_DFTB

   implicit none
   integer, intent(in)    :: M_in
   integer, intent(in)    :: M
   real*8 , intent(in)    :: Smat(M,M)
   real*8 , intent(inout) :: X_min(M,M)
   real*8 , intent(out)   :: Xtrans(M_in,M_in), Ytrans(M_in,M_in), &
                             Xmm(M_in,M_in),eigenvalues(M)
   real*8 , intent(inout) :: Xmat(M_in,M_in), Ymat(M_in,M_in)
   real*8                 :: Y_min(M,M)

   integer :: icount, jcount, LWORK, info
   real*8, allocatable :: WORK(:)

   ! Diagonalization of S matrix, both with ESSL or LAPACK.
   ! The S matrix is stored in RMM(M13, M13+1, ..., M13+MM).
   do icount = 1, M
   do jcount = 1, M
      X_min(icount, jcount) = Smat(icount, jcount)
   enddo
   enddo

   if (allocated(WORK)) deallocate(WORK)
   allocate(WORK(1))
   call dsyev('V', 'L', M, X_min, M, eigenvalues, WORK, -1, info)

   LWORK = int(WORK(1))
   deallocate(WORK)
   allocate(WORK(LWORK))
   call dsyev('V', 'L', M, X_min, M, eigenvalues, WORK, LWORK, info)


   ! Here we obtain the transformation matrices X and Y for converting from the
   ! atomic orbital basis to the molecular orbital basis (truncated during
   ! linear dependency elimination). S is the overlap matrix, s is the diagonal
   ! eigenvalue matrix of S and U is the eigenvector matrix of S:
   ! X = U s^(-1/2)
   ! Matrix X's dimension is M*3M. In the first M*M terms it contains the
   ! transformation matrices and in the other M*2M terms it contains auxiliar
   ! matrices.
   do jcount = 1, M
      if (eigenvalues(jcount).lt.1.0D-06) then
         write(*,*) 'WARNING - TD: Linear dependency detected in S matrix.'
         do icount = 1, M
            X_min(icount, jcount) = 0.0D0
            Y_min(icount, jcount) = 0.0D0
         enddo
      else
         do icount = 1, M
            Y_min(icount,jcount) = X_min(icount,jcount) *sqrt(eigenvalues(jcount))
            X_min(icount,jcount) = X_min(icount,jcount) /sqrt(eigenvalues(jcount))
         enddo
      endif
   enddo

!------------------------------------------------------------------------------!
!DFTB: Xmat and Ymat are adapted for DFTB
      if (dftb_calc) then
         call getXY_DFTB(M,X_min,Y_min,Xmat,Ymat)
      else

         xmat=x_min(1:M,1:M)
         ymat=y_min

      end if

!------------------------------------------------------------------------------!
   ! Stores transformation matrix Xmm and transposed matrices.
   do icount = 1, M_in
   do jcount = 1, M_in
      Xmm(icount, jcount)    = Xmat(icount, jcount)
      Xtrans(jcount, icount) = Xmat(icount, jcount)
      Ytrans(jcount, icount) = Ymat(icount, jcount)
   enddo
   enddo

   return
end subroutine td_overlap_diag

subroutine td_bc_rho(M, rho, Ymat, Ytrans)
   implicit none
   integer, intent(in) :: M
   real*8 , intent(in) :: Ymat(M,M), Ytrans(M,M)
#ifdef TD_SIMPLE
   complex*8 , intent(inout) :: rho(M,M)
#else
   complex*16, intent(inout) :: rho(M,M)
#endif

   rho = matmul(Ytrans, rho)
   rho = matmul(rho   , Ymat)

   return
end subroutine td_bc_rho

subroutine td_coulomb_precalc(igpu, MEMO)
   use faint_cpu77, only: int3mem
   implicit none
   integer, intent(in)    :: igpu
   logical, intent(inout) :: MEMO

   if (igpu.gt.2) then
      call g2g_timer_start('Coulomb - precalc')
      call aint_coulomb_init()
      if (igpu.eq.5) MEMO = .false.
   endif
   if (MEMO) then
      call int3mem()
   endif
   call g2g_timer_stop('Coulomb - precalc')

   return
end subroutine td_coulomb_precalc

! Subroutines used in propagation.
subroutine td_get_time(t, tdstep, istep, propagator, is_lpfrg)
   implicit none
   integer, intent(in)    :: istep, propagator
   real*8 , intent(in)    :: tdstep
   logical, intent(in)    :: is_lpfrg
   real*8 , intent(inout) :: t

   select case (propagator)
      case (1)
         t = (istep-1) * tdstep
      case (2)
         if (is_lpfrg) then
            t = (istep-1) * tdstep * 0.1
         else
            t = 20 * tdstep
            t = t + (istep-200)*tdstep
         endif
      case default
   end select

   t = t * 0.024190D0
   write(*,*) 'TD - Time (fs)  =', t
   return
end subroutine td_get_time

subroutine td_check_prop(is_lpfrg, propagator, istep, lpfrg_steps, tdrestart)
   implicit none
   integer, intent(in)  :: propagator, istep, lpfrg_steps
   logical, intent(in)  :: tdrestart
   logical, intent(out) :: is_lpfrg

   is_lpfrg = ((propagator.eq.1) .or. (((propagator.eq.2) .and. &
              (istep.lt.lpfrg_steps)) .and. (.not.tdrestart)))
   if ( (is_lpfrg).and.(istep.eq.1) ) then
      write(*,*) 'TD - Starting Verlet Propagation'
   endif
   if ( (.not.(is_lpfrg)).and.(((istep-1).eq.lpfrg_steps)) ) then
      write(*,*) 'TD - Starting Magnus Propagation'
   endif
   return
end subroutine td_check_prop

subroutine td_calc_energy(E, E1, E2, En, Ex, Es, MM, RMM, RMM11, is_lpfrg, &
                          transport_calc, sol, time)
   use faint_cpu77, only: int3lu
   use field_subs , only: field_calc
   implicit none
   integer, intent(in)    :: MM
   logical, intent(in)    :: is_lpfrg, transport_calc, sol
   real*8 , intent(in)    :: time
   real*8 , intent(inout) :: E, E1, E2, En, Ex, Es, RMM(MM), RMM11(MM)
   integer :: icount

   E1 = 0.0D0; E = 0.0D0
   if (is_lpfrg) then
      call g2g_timer_sum_start("TD - Coulomb")
      call int3lu(E2)
      call g2g_timer_sum_pause("TD - Coulomb")
      call g2g_timer_sum_start("TD - Exc")
      call g2g_solve_groups(0,Ex,0)
      call g2g_timer_sum_pause("TD - Exc")
   endif

   ! ELECTRIC FIELD CASE - Perturbation type: Gaussian.
   if (.not.transport_calc) call field_calc(E1, time)

   ! Add 1e contributions to E1.
   do icount = 1, MM
      E1 = E1 + RMM(icount)*RMM11(icount)
   enddo
   E = E1 + E2 + En + Ex
   if (sol) E = E + Es

   return
end subroutine td_calc_energy

subroutine td_bc_fock(M_in, M, MM, RMM5, fock, Xtrans, Xmm, natom, nshell,ncont, istep)

   use dftb_data, only:dftb_calc,MTB
   use dftb_subs, only:chimeraDFTB_evol
   implicit none
   integer, intent(in)    :: M, MM, M_in
   real*8 , intent(inout) :: RMM5(MM), Xtrans(M_in,M_in), Xmm(M_in,M_in)
   real*8, intent(out)    :: fock(M_in, M_in)
   real*8                 :: Xtemp(M_in,M_in)
   real*8                 :: fock_0(M,M)
   integer, intent(in)  :: natom
   integer, intent(in)  :: ncont(M)
   integer, intent(in)  :: istep
   integer, intent(in)  :: nshell (0:4)

   call g2g_timer_start('fock')
   call spunpack('L', M, RMM5, fock_0)

   if (dftb_calc) then
      call chimeraDFTB_evol(M,fock_0, fock, natom, nshell,ncont, istep)
   else

      fock=fock_0

   end if

   Xtemp = matmul(Xtrans, fock)
   fock  = matmul(Xtemp , Xmm)
   call sprepack('L', M, RMM5, fock(MTB+1:MTB+M,MTB+1:MTB+M))
   call g2g_timer_stop('fock')

   return
end subroutine td_bc_fock

subroutine td_verlet(M, M_in, fock, rhold, rho, rhonew, istep, Im, dt_lpfrg,         &
                     transport_calc, natom, Nuc, Iz, overlap, sqsm, rho_aux,         &
                     Xmat, Xtrans)
   use mathsubs      , only : commutator
   use transport_subs, only : transport_propagate
   use dftb_data     , only : dftb_calc, rhold_AOTB, rhonew_AOTB
   implicit none
   integer   , intent(in)    :: M, M_in, istep, natom, Nuc(M), Iz(natom)
   real*8    , intent(in)    :: dt_lpfrg
   real*8    , intent(in)    :: Xmat(M_in, M_in), Xtrans(M_in, M_in)
   logical   , intent(in)    :: transport_calc
   real*8    , intent(inout) :: overlap(:,:), sqsm(M,M), fock(M_in,M_in)
#ifdef TD_SIMPLE
   complex*8 , intent(inout) :: Im, rho_aux(M,M), rhold(M_in,M_in),   &
                                rhonew(M_in,M_in), rho(M_in,M_in)
#else
   complex*16, intent(inout) :: Im, rho_aux(M_in,M_in), rhold(M_in,M_in),   &
                                rhonew(M_in,M_in), rho(M_in,M_in)
#endif
   integer :: icount, jcount

   if (istep.eq.1) then
      call g2g_timer_start('conmutc')
      rhold = commutator(fock,rho)
      rhold = rho + dt_lpfrg*(Im*rhold)
      call g2g_timer_stop('conmutc')
   endif

   if ((transport_calc) .and. (istep.ge.3))then
      call transport_propagate(M, natom, Nuc, Iz, 1, istep, overlap, sqsm, &
                               rho_aux)
   endif

   call g2g_timer_start('commutator')
   rhonew = commutator(fock, rho)
   rhonew = rhold - dt_lpfrg*(Im*rhonew)
   call  g2g_timer_stop('commutator')

   !Transport: Add the driving term to the propagation.
   if ((istep.ge.3) .and. (transport_calc)) then
      write(*,*) 'Transport: Adding driving term to the density.'
      rhonew = rhonew - rho_aux
   endif

   ! Density update (rhold-->rho, rho-->rhonew)
   do icount = 1, M_in
   do jcount = 1, M_in
      rhold(icount, jcount) = rho(icount, jcount)
      rho(icount, jcount)   = rhonew(icount, jcount)
   enddo
   enddo

!DFTB: rhonew and rhold in AO is store for charge calculations of DFTB
   if (dftb_calc) then

      rhold_AOTB=matmul(xmat,rhold)
      rhold_AOTB=matmul(rhold_AOTB,xtrans)

      rhonew_AOTB=matmul(xmat,rhonew)
      rhonew_AOTB=matmul(rhonew_AOTB,xtrans)

   end if

   return
end subroutine td_verlet

subroutine td_magnus(M, fock, F1a, F1b, rho, rhonew, factorial, NBCH,  &
                     dt_magnus, natom, transport_calc, Nuc, Iz, istep, &
                     overlap, sqsm, rho_aux, Xmat, Xtrans, time, M_in, &
                     nshell, ncont)
   use transport_subs, only: transport_propagate
   use dftb_data,      only:dftb_calc,MTB, rhold_AOTB, rhonew_AOTB
   use dftb_subs,      only:chimeraDFTB_evol
   implicit none
   logical  , intent(in)     :: transport_calc
   integer  , intent(in)     :: M, NBCH, istep, natom, Nuc(M), Iz(natom), M_in
   real*8   , intent(in)     :: dt_magnus, factorial(NBCH), Xtrans(M,M), time
   real*8   , intent(inout)  :: fock(M_in,M_in), F1a(M_in,M_in), &
                                F1b(M_in,M_in), overlap(:,:), sqsm(M,M), &
                                Xmat(M,M)
   integer, intent(in)       :: nshell (0:4)
   integer, intent(in)       :: ncont(M)
#ifdef TD_SIMPLE
   complex*8 , intent(inout) :: rho_aux(M_in,M_in), rhonew(M_in,M_in), rho(M_in,M_in)
#else
   complex*16, intent(inout) :: rho_aux(M_in,M_in), rhonew(M_in,M_in), rho(M_in,M_in)
#endif
   real*8                    :: fock_aux(M_in,M_in)

   if (transport_calc) then
      call transport_propagate(M, natom, Nuc, Iz, 2, istep, overlap, sqsm,     &
                               rho_aux)
   endif

! DFTB: this if is temporary, it is to conserve the atomic part of dftb in
! predictor.
   if (dftb_calc) then
      call chimeraDFTB_evol(M,fock(MTB+1:MTB+M,MTB+1:MTB+M), fock_aux, natom,      &
                        nshell,ncont, istep)
      fock=fock_aux

!DFTB: rhold in AO is store for charge calculations of DFTB
      rhold_AOTB=matmul(xmat,rho)
      rhold_AOTB=matmul(rhold_AOTB,xtrans)
   end if

   call g2g_timer_start('predictor')
   call predictor(F1a, F1b, fock, rho, factorial, Xmat, Xtrans, dt_magnus, &
                  time, M_in, MTB)
   call g2g_timer_stop('predictor')
   call g2g_timer_start('magnus')
   call magnus(fock, rho, rhonew, M_in, NBCH, dt_magnus, factorial)
   call g2g_timer_stop('magnus')

   !Transport: Add the driving term to the propagation.
   if (transport_calc) then
      write(*,*) 'Transport: Adding driving term to the density.'
      rhonew = rhonew - rho_aux
   endif
!DFTB: rhonew in AO is store for charge calculations of DFTB
   if (dftb_calc) then
      rhonew_AOTB=matmul(xmat,rhonew)
      rhonew_AOTB=matmul(rhonew_AOTB,xtrans)
   end if

   ! Density update and Fock storage.
   F1a = F1b
   F1b = fock
   rho = rhonew
   return
end subroutine td_magnus

subroutine td_dipole(t, tdstep, Fx, Fy, Fz, istep, propagator, is_lpfrg, uid)
   implicit none
   integer, intent(in)    :: istep, propagator, uid
   logical, intent(in)    :: is_lpfrg
   real*8 , intent(in)    :: Fx, Fy, Fz, t, tdstep
   real*8 :: dipxyz(3)

   if(istep.eq.1) then
      call write_dipole_td_header(tdstep, Fx, Fy, Fz, uid)
   endif
   if ((propagator.gt.1).and.(is_lpfrg)) then
      if (mod ((istep-1),10) == 0) then
         call g2g_timer_start('DIPOLE_TD')
         call dip(dipxyz)
         call g2g_timer_stop('DIPOLE_TD')
         call write_dipole_td(dipxyz, t, uid)
      endif
   else
      call g2g_timer_start('DIPOLE_TD')
      call dip(dipxyz)
      call g2g_timer_stop('DIPOLE_TD')
      call write_dipole_td(dipxyz, t, uid)
   endif

   return
end subroutine td_dipole

! CUBLAS-dependent subroutines.
#ifdef CUBLAS
subroutine td_allocate_cublas(M, sizeof_real, devPtrX, devPtrY)
   implicit none
   external CUBLAS_INIT
   integer  , intent(in)    :: M, sizeof_real
   integer*8, intent(inout) :: devPtrX, devPtrY
   integer :: stat
   stat = 0

   write(*,*) 'TD: Using CUBLAS.'
   call CUBLAS_INIT()
   return
end subroutine td_allocate_cublas

subroutine td_finalise_cublas(devPtrX, devPtrY, devPtrXc)
   implicit none
   external CUBLAS_SHUTDOWN, CUBLAS_FREE
   integer*8 :: devPtrX, devPtrXc, devPtrY
   call CUBLAS_FREE(devPtrX)
   call CUBLAS_FREE(devPtrXc)
   call CUBLAS_FREE(devPtrY)
   call CUBLAS_SHUTDOWN()
end subroutine td_finalise_cublas

subroutine td_bc_rho_cu(M, sizeof_real, sizeof_complex, devPtrX, devPtrY,      &
                        devPtrXc, Xmat, Ymat, rho, rho_aux)
   use cublasmath, only : basechange_cublas
   implicit none
   external CUBLAS_SHUTDOWN, CUBLAS_ALLOC, CUBLAS_SET_MATRIX
   integer  , intent(in) :: M, sizeof_real, sizeof_complex
   integer*8, intent(inout) :: devPtrX, devPtrY, devPtrXc
   real*8   , intent(inout) :: Ymat(M,M), Xmat(M,M)
#ifdef TD_SIMPLE
   complex*8 , intent(inout) :: rho(M,M), rho_aux(M,M)
#else
   complex*16, intent(inout) :: rho(M,M), rho_aux(M,M)
#endif
   integer :: icount, jcount, stat

   do icount = 1, M
   do jcount = 1, M
      rho_aux(icount, jcount) = cmplx(Xmat(icount, jcount), 0.0D0)
   enddo
   enddo

   call CUBLAS_ALLOC(M*M, sizeof_real   , devPtrX)
   call CUBLAS_ALLOC(M*M, sizeof_complex, devPtrXc)
   call CUBLAS_ALLOC(M*M, sizeof_complex, devPtrY)
   if (stat.NE.0) then
      write(*,*) "ERROR - TD: X and/or Y CUBLAS memory allocation failed."
      call CUBLAS_SHUTDOWN()
      stop
   endif

   call CUBLAS_SET_MATRIX(M, M, sizeof_complex, rho_aux, M, devPtrXc, M)
   call CUBLAS_SET_MATRIX(M, M, sizeof_real   , Xmat   , M, devPtrX , M)
   do icount = 1, M
   do jcount = 1, M
      rho_aux(icount, jcount) = cmplx(Ymat(icount, jcount), 0.0D0)
   enddo
   enddo


   call CUBLAS_SET_MATRIX(M, M, sizeof_complex, rho_aux, M, devPtrY, M)
   if (stat.NE.0) then
      write(*,*) "ERROR - TD : X and/or Y CUBLAS setting failed."
      call CUBLAS_SHUTDOWN
      stop
   endif
   rho_aux = 0.0D0
   call g2g_timer_start('complex_rho_ao_to_on-cu')
   rho_aux = basechange_cublas(M, rho, devPtrY, 'dir')
   rho     = rho_aux
   call g2g_timer_stop('complex_rho_ao_to_on-cu')

   return
end subroutine td_bc_rho_cu

subroutine td_bc_fock_cu(M_in,M, MM, RMM5, fock, devPtrX, natom, nshell, ncont, istep)
   use dftb_data, only:dftb_calc, MTB
   use dftb_subs, only:chimeraDFTB_evol
   use cublasmath, only: basechange_cublas
   implicit none
   integer  , intent(in)    :: M, MM, M_in
   integer*8, intent(in)    :: devPtrX
   real*8   , intent(inout) :: RMM5(MM)
   real*8   , intent(out)   :: fock(M_in,M_in)
   real*8 :: fock_0(M,M)
   real*8 :: Xtemp(M_in,M_in)
   integer, intent(in)  :: natom
   integer, intent(in)  :: ncont(M)
   integer, intent(in)  :: istep
   integer, intent(in)  :: nshell (0:4)

   call g2g_timer_start('fock')
   call spunpack('L', M, RMM5, fock_0)

   if (dftb_calc) then
      call chimeraDFTB_evol(M,fock_0, fock, natom, nshell,ncont, istep)
   else
      fock=fock_0
   end if

   Xtemp = basechange_cublas(M_in, fock, devPtrX, 'dir')
   fock  = Xtemp
   call sprepack('L', M, RMM5, fock(MTB+1:MTB+M,MTB+1:MTB+M))
   call g2g_timer_stop('fock')

   return
end subroutine td_bc_fock_cu

subroutine td_verlet_cu(M, M_in, fock, rhold, rho, rhonew, istep, Im, dt_lpfrg,      &
                     transport_calc, natom, Nuc, Iz, overlap, sqsm, rho_aux,   &
                     devPtrY, devPtrXc)
   use cublasmath    , only : commutator_cublas,basechange_cublas
   use transport_subs, only : transport_propagate_cu
   use dftb_data     , only : dftb_calc, rhold_AOTB, rhonew_AOTB
   implicit none
   integer   , intent(in)    :: M,M_in, istep, natom, Nuc(M), Iz(natom)
   integer*8 , intent(in)    :: devPtrY, devPtrXc
   real*8    , intent(in)    :: dt_lpfrg
   logical   , intent(in)    :: transport_calc
   real*8    , intent(inout) :: overlap(:,:), sqsm(M,M), fock(M_in,M_in)
#ifdef TD_SIMPLE
   complex*8 , intent(inout) :: Im, rho_aux(M,M), rhold(M_in,M_in), rhonew(M_in,M_in), &
                                rho(M_in,M_in)
#else
   complex*16, intent(inout) :: Im, rho_aux(M_in,M_in), rhold(M_in,M_in), rhonew(M_in,M_in), &
                                rho(M_in,M_in)
#endif
   integer :: icount, jcount

   if(istep.eq.1) then
      call g2g_timer_start('cuconmut')
      rhold = commutator_cublas(fock, rho)
      call g2g_timer_stop('cuconmut')
      rhold = rho + dt_lpfrg*(Im*rhold)
   endif

   if ((transport_calc) .and. (istep.ge.3)) then
      call transport_propagate_cu(M, natom, Nuc, Iz, 1, istep, overlap, sqsm, &
                                  rho_aux, devPtrY)
   endif

   call g2g_timer_start('commutator')
   rhonew = commutator_cublas(fock, rho)
   rhonew = rhold - dt_lpfrg*(Im*rhonew)
   call  g2g_timer_stop('commutator')

   !Transport: Add the driving term to the propagation.
   if ((istep.ge.3) .and. (transport_calc)) then
      write(*,*) 'Transport: Adding driving term to the density.'
      rhonew = rhonew - rho_aux
   endif

   ! Density update (rhold-->rho, rho-->rhonew)
   do icount = 1, M_in
   do jcount = 1, M_in
      rhold(icount, jcount) = rho(icount,jcount)
      rho(icount,jcount)    = rhonew(icount,jcount)
   enddo
   enddo

!DFTB: rhonew and rhold in AO is store for charge calculations of DFTB
   if (dftb_calc) then
      rhold_AOTB=basechange_cublas(M_in,rhold,devPtrXc,'inv')
      rhonew_AOTB=basechange_cublas(M_in,rhonew,devPtrXc,'inv')
   end if


end subroutine td_verlet_cu

subroutine td_magnus_cu(M, fock, F1a, F1b, rho, rhonew, devPtrX, devPtrXc,     &
                        factorial, NBCH, dt_magnus, natom, transport_calc, Nuc,&
                        Iz, istep, overlap, sqsm, rho_aux, devPtrY, time, M_in,&
                        nshell,ncont)
   use cublasmath    , only: cupredictor, cumagnusfac, basechange_cublas
   use transport_subs, only: transport_propagate_cu
   use dftb_data,      only:dftb_calc,MTB, rhold_AOTB, rhonew_AOTB
   use dftb_subs,      only:chimeraDFTB_evol
   implicit none
   logical  , intent(in)     :: transport_calc
   integer  , intent(in)     :: M, NBCH, istep, natom, Nuc(M), Iz(natom)
   integer  , intent(in)     :: M_in
   integer*8, intent(in)     :: devPtrX, devPtrXc, devPtrY
   real*8   , intent(in)     :: dt_magnus, factorial(NBCH), time
   integer  , intent(in)     :: nshell (0:4), ncont(M)
   real*8   , intent(inout)  :: fock(M_in,M_in), F1a(M_in,M_in), &
                                F1b(M_in,M_in), overlap(:,:), sqsm(M,M)
#ifdef TD_SIMPLE
   complex*8 , intent(inout) :: rho_aux(M_in,M_in), rhonew(M_in,M_in), rho(M_in,M_in)
#else
   complex*16, intent(inout) :: rho_aux(M_in,M_in), rhonew(M_in,M_in), rho(M_in,M_in)
#endif
   real*8   :: fock_aux(M_in, M_in)

   if (transport_calc) then
      call transport_propagate_cu(M, natom, Nuc, Iz, 2, istep, overlap, sqsm, &
                                  rho_aux, devPtrY)
   endif

! DFTB: this if is temporary, it is to conserve the atomic part of dftb in
! predictor.
   if (dftb_calc) then
      call chimeraDFTB_evol(M,fock(MTB+1:MTB+M,MTB+1:MTB+M), fock_aux, natom,      &
                        nshell,ncont, istep)
      fock=fock_aux
      !DFTB: rhold in AO is store for charge calculations of DFTB
      rhold_AOTB=basechange_cublas(M_in,rho,devPtrXc,'inv')
   end if

   call g2g_timer_start('cupredictor')
   call cupredictor(F1a, F1b, fock, rho, devPtrX, factorial, devPtrXc, &
                    dt_magnus, time, M_in, MTB)
   call g2g_timer_stop('cupredictor')
   call g2g_timer_start('cumagnus')
   call cumagnusfac(fock, rho, rhonew, M_in, NBCH, dt_magnus, factorial)
   call g2g_timer_stop('cumagnus')

   !Transport: Add the driving term to the propagation.
   if (transport_calc) then
      write(*,*) 'Transport: Adding driving term to the density.'
      rhonew = rhonew - rho_aux
   endif

!DFTB: rhonew in AO is store for charge calculations of DFTB
   if (dftb_calc) then
      rhonew_AOTB=basechange_cublas(M_in,rhonew,devPtrXc,'inv')
   end if

   ! Density update and Fock storage.
   F1a = F1b
   F1b = fock
   rho = rhonew
   return
end subroutine td_magnus_cu

#endif

end module time_dependent
