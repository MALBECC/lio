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

subroutine TD(fock_aop, rho_aop, fock_bop, rho_bop)
   use garcha_mod    , only: M, Md, NBCH, propagator, RMM, NCO, Iz, igrid2, r, &
                             Nuc, nsol, pc, X, Smat, MEMO, sol, natom, sqsm,   &
                             Nunp, ntatom, ncont, nshell, a, c, d, NORM,OPEN, &
                             rhoalpha, rhobeta
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
                             read_td_restart_verlet , read_td_restart_magnus,  &
                             write_energies
   use fileio_data   , only: verbose
   use typedef_operator, only: operator

#ifdef CUBLAS
   use cublasmath    , only: basechange_cublas
#endif

   implicit none
!carlos: Operator inserted for TD
   type(operator), intent(inout)           :: rho_aop, fock_aop
   type(operator), intent(inout), optional :: rho_bop, fock_bop

   real*8  :: E, En, E1, E2, E1s, Es, Ens = 0.0D0, Ex, t, dt_magnus, dt_lpfrg
   integer :: MM, MMd, M2, M3 ,M5, M13, M15, M11, LWORK, igpu, info, istep,    &
              icount,jcount
   integer :: lpfrg_steps = 200, chkpntF1a = 185, chkpntF1b = 195
   logical :: is_lpfrg
   character(len=20) :: restart_filename

   real*8 , allocatable, dimension(:)   :: factorial, WORK
   real*8 , allocatable, dimension(:,:) :: Xmm, Xmat,Xtrans, Ytrans, overlap,  &
                                           Ymat
!carlos: the next variables have 3 dimensions, the 3th one is asociated with the
!        spin number. This one will have the value of 1 for Close Shell and 2
!        for Open shell. Spin=1 will always be reffered to alpha and spin=2,
!        beta.

   real*8 , allocatable, dimension(:,:,:) :: fock, F1a, F1b
   real*8 , allocatable, dimension(:,:,:) :: fock_0

! Precision options.
#ifdef TD_SIMPLE
   complex*8  :: Im = (0.0E0,2.0E0)
   complex*8 , allocatable, dimension(:,:,:) :: rho, rho_aux, rhonew, rhold
   complex*8, allocatable, dimension(:,:,:) :: rho_0
#else
   complex*16 :: Im = (0.0D0,2.0D0)
   complex*16, allocatable, dimension(:,:,:) :: rho, rho_aux, rhonew, rhold
   complex*16, allocatable, dimension(:,:,:) :: rho_0
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
!carlos: Open Shell variables

   integer :: NCOa, NCOb
   integer :: dim3
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
   M3  = 1+MM
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
!carlos: dim3 is the 3th dimension of all matrix involve in open shell
!        calculation

   NCOa = NCO
   dim3 = 1
   if (open) then
     NCOb = NCO + Nunp
     dim3 = 2
   end if

!------------------------------------------------------------------------------!

   ! Checks and performs allocations.
   call td_allocate_all(M_in, M, dim3, NBCH, propagator, F1a, F1b, fock, rho,  &
                        rho_aux, rhold, rhonew, rho_0, fock_0, sqsm, Xmm, Xmat,&
                        Xtrans, Ymat, Ytrans, factorial)

   ! Initialises propagator-related parameters and other variables.
   call td_initialise(propagator, tdstep, NBCH, dt_lpfrg, dt_magnus, factorial,&
                      NCO, Nunp, natom, Iz)
!DFTB:TD restart is not working with DFTB
   ! TD restart reading.
   if (tdrestart) then
      if (OPEN) then
         restart_filename = 'td_a_in.restart'
         if (propagator.eq.2) then
            call read_td_restart_magnus(rho(:,:,1), F1a(:,:,1), F1b(:,:,1),    &
                                        M_in, restart_filename)
         else
            call read_td_restart_verlet(rho(:,:,1), M_in, restart_filename)
         endif
         restart_filename = 'td_b_in.restart'
         if (propagator.eq.2) then
            call read_td_restart_magnus(rho(:,:,2), F1a(:,:,2), F1b(:,:,2),    &
                                        M_in, restart_filename)
         else
            call read_td_restart_verlet(rho(:,:,2), M_in, restart_filename)
         endif
         call sprepack_ctr('L', M, rhoalpha, rho(MTB+1:MTB+M,MTB+1:MTB+M,1))
         call sprepack_ctr('L', M, rhobeta, rho(MTB+1:MTB+M,MTB+1:MTB+M,2))
         RMM(1:MM) = rhoalpha + rhobeta
      else
         restart_filename = 'td_in.restart'
         if (propagator.eq.2) then
            call read_td_restart_magnus(rho(:,:,1), F1a(:,:,1), F1b(:,:,1),    &
                                        M_in, restart_filename)
         else
            call read_td_restart_verlet(rho(:,:,1), M_in, restart_filename)
         endif
         call sprepack_ctr('L', M, RMM, rho(MTB+1:MTB+M,MTB+1:MTB+M,1))
      end if
   else
      ! Read the density matrix stored in RMM(1,2,3,...,MM) into rho matrix.
      !For Open Shell, this is read from rhoalpha and rhobeta
      if (OPEN) then
         call spunpack_rtc('L', M, rhoalpha, rho_0(:,:,1))
         call spunpack_rtc('L', M, rhobeta, rho_0(:,:,2))
      else
         call spunpack_rtc('L', M, RMM, rho_0(:,:,1))
       end if
   endif

!------------------------------------------------------------------------------!
! DFTB: Initialize DFTB variables for TD
   if (dftb_calc) then
      call dftb_td_init (M , rho, rho_0, overlap, RMM(M5),dim3)
   else
      rho=rho_0
   end if
!carlos: storing rho AO data in Operator

   call rho_aop%Sets_dataC_AO(rho(:,:,1))

   if (OPEN) call rho_bop%Sets_dataC_AO(rho(:,:,2))

!------------------------------------------------------------------------------!

   ! Proper TD calculation start.
   if (verbose .gt. 2) write(*,'(A)') 'Starting TD calculation.'
   ! Create integration grid for XC, assigning points to groups (spheres/cubes)
   ! and significant functions to groups, also calculating point weights.
   if (field) call field_setup_old(pert_time, 1, fx, fy, fz)
   call td_integration_setup(igrid2, igpu)
   call td_integral_1e(E1, En, E1s, Ens, MM, igpu, nsol, RMM, RMM(M11), r, pc, &
                       ntatom,natom,Smat,Nuc,a,c,d,Iz,ncont,NORM,M,Md,nshell)
   ! Initialises transport if required.
   if (transport_calc) call transport_init(M, dim3, natom, Nuc, RMM(M5),       &
                                           overlap, rho,OPEN)

   ! Diagonalizes Smat, stores transformation matrix Xmm and calculates the
   ! transposed matrices Xtrans and Ytrans
   call td_overlap_diag(M_in, M, Smat, RMM(M13), X, Xmat ,Xtrans, Ymat, Ytrans,&
                        Xmm)
   ! Here rho_aux is used as an auxiliar matrix for CUBLAS. Then we need to
   ! restore its original value. Rho is then transformed to the orthonormal
   ! basis with either matmul intrinsic or CUBLAS.

!carlos: td_allocate_cublas now initialize CUBLAS and allocates the memory to
!        Xmat and Ymat necessary for base changes
#ifdef CUBLAS
   call td_allocate_cublas(M_in, sizeof_real, sizeof_complex, devPtrX, devPtrY,&
                          devPtrXc, Xmat, Ymat, rho_aux(:,:,1))
   if (verbose .gt. 3) write(*,'(A)') "  TD - Using CUBLAS"
#endif

!carlos: now the base changes is inside the operatros

#ifdef CUBLAS
   call rho_aop%BChange_AOtoON(devPtrY, M_in, 'c')
   if (OPEN) call rho_bop%BChange_AOtoON(devPtrY, M_in, 'c')
#else
   call rho_aop%BChange_AOtoON(Ymat, M_in , 'c')
   if(OPEN) call rho_bop%BChange_AOtoON(Ymat, M_in , 'c')
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
      call td_check_prop(is_lpfrg, propagator, istep, lpfrg_steps, tdrestart, &
                         verbose)
      call td_get_time(t, tdstep, istep, propagator, is_lpfrg)

      call g2g_timer_sum_start("TD - TD Step Energy")

      call td_calc_energy(E, E1, E2, En, Ex, Es, MM, RMM, RMM(M11), is_lpfrg,  &
                          transport_calc, sol, t/0.024190D0)
      call g2g_timer_sum_pause("TD - TD Step Energy")
      if (verbose .gt. 2) write(*,'(A,I6,A,F12.6,A,F12.6,A)') "  TD Step: ", &
                                istep, " - Time: ", t, " fs - Energy : ", E, &
                                " A.U."

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
         call td_bc_fock_cu(M_in,M, MM, RMM(M5), fock_aop, devPtrX,natom,      &
                            nshell,ncont, istep)

         if (OPEN) then

            call td_bc_fock_cu(M_in,M, MM, RMM(M3), fock_bop, devPtrX, natom,  &
                               nshell,ncont, istep)
            call td_verlet_cu(M, M_in, dim3, OPEN, fock_aop, rhold, rho_aop,   &
                              rhonew, istep, Im, dt_lpfrg, transport_calc,     &
                              natom, Nuc, Iz, overlap, sqsm, devPtrY,devPtrXc, &
                              fock_bop, rho_bop)
         else

            call td_verlet_cu(M, M_in, dim3, OPEN, fock_aop, rhold, rho_aop,   &
                              rhonew, istep, Im, dt_lpfrg, transport_calc,     &
                              natom, Nuc, Iz, overlap, sqsm, devPtrY,devPtrXc)
         end if

         if (propagator.eq.2) then
            call fock_aop%Gets_data_ON(fock(:,:,1))
            if (OPEN) call fock_bop%Gets_data_ON(fock(:,:,2))
            if (istep.eq.chkpntF1a) F1a = fock
            if (istep.eq.chkpntF1b) F1b = fock
         endif
      else

         if(OPEN) then
            call td_magnus_cu(M, dim3, OPEN,fock_aop, F1a, F1b, rho_aop,       &
                              rhonew,devPtrX, devPtrXc, factorial, NBCH,       &
                              dt_magnus, natom, transport_calc, Nuc, Iz,       &
                              istep, overlap,sqsm, devPtrY,                    &
                              t/0.024190D0, M_in, nshell, ncont, fock_bop,     &
                              rho_bop)
         else
            call td_magnus_cu(M, dim3, OPEN, fock_aop, F1a, F1b, rho_aop,      &
                              rhonew, devPtrX, devPtrXc, factorial, NBCH,      &
                              dt_magnus, natom, transport_calc, Nuc, Iz,       &
                              istep, overlap, sqsm, devPtrY,                   &
                              t/0.024190D0, M_in, nshell, ncont)
         endif
      end if

      call g2g_timer_start('complex_rho_on_to_ao-cu')
      call rho_aop%BChange_ONtoAO(devPtrXc, M_in, 'c')
      if (OPEN) call rho_bop%BChange_ONtoAO(devPtrXc, M_in, 'c')
      call g2g_timer_stop('complex_rho_on_to_ao-cu')

#else
      if (is_lpfrg) then
         call td_bc_fock(M_in, M, MM, RMM(M5), fock_aop,Xmat, natom, nshell,    &
                         ncont, istep)
         if (OPEN) then
            call td_bc_fock(M_in, M, MM, RMM(M3), fock_bop,Xmat, natom, nshell, &
                            ncont, istep)

            call td_verlet(M, M_in, dim3, OPEN, fock_aop, rhold, rho_aop,      &
                           rhonew, istep, Im, dt_lpfrg, transport_calc, natom, &
                           Nuc, Iz, overlap, sqsm, Xmat, Ymat, Xtrans,         &
                           fock_bop, rho_bop)
         else
            call td_verlet(M, M_in, dim3, OPEN, fock_aop, rhold, rho_aop,      &
                           rhonew, istep, Im, dt_lpfrg, transport_calc, natom, &
                           Nuc, Iz, overlap, sqsm, Xmat, Ymat, Xtrans)
         end if

         if (propagator.eq.2) then
            call fock_aop%Gets_data_ON(fock(:,:,1))
            if (OPEN) call fock_bop%Gets_data_ON(fock(:,:,2))
            if (istep.eq.chkpntF1a) F1a = fock
            if (istep.eq.chkpntF1b) F1b = fock
         endif
      else

         if(OPEN) then
            call td_magnus(M, dim3, OPEN, fock_aop, F1a, F1b, rho_aop, rhonew, &
                           factorial,NBCH, dt_magnus, natom, transport_calc,   &
                           Nuc, Iz, istep, overlap, sqsm, Xmat, Ymat, Xtrans,  &
                           t/0.024190D0, M_in, nshell, ncont, fock_bop, rho_bop)
         else
            call td_magnus(M, dim3, OPEN, fock_aop, F1a, F1b, rho_aop, rhonew, &
                           factorial,NBCH, dt_magnus, natom, transport_calc,   &
                           Nuc, Iz, istep, overlap, sqsm, Xmat, Ymat, Xtrans,  &
                           t/0.024190D0, M_in, nshell, ncont)
         end if
      endif

      call g2g_timer_start('complex_rho_on_to_ao')
      call rho_aop%BChange_ONtoAO(Xtrans, M_in, 'c')
      if (OPEN) call rho_bop%BChange_ONtoAO(Xtrans, M_in, 'c')
      call g2g_timer_stop('complex_rho_on_to_ao')
#endif
      call g2g_timer_sum_pause("TD - Propagation")

      ! The real part of the density matrix in the atomic orbital basis is
      ! copied in RMM(1,2,3,...,MM) to compute the corresponding Fock matrix.
      if (OPEN) then
         call rho_aop%Gets_dataC_AO(rho_aux(:,:,1))
         call rho_bop%Gets_dataC_AO(rho_aux(:,:,2))
         call sprepack_ctr('L', M, rhoalpha, rho_aux(MTB+1:MTB+M,MTB+1:MTB+M,1))
         call sprepack_ctr('L', M, rhobeta, rho_aux(MTB+1:MTB+M,MTB+1:MTB+M,2))

         RMM(1:MM) = rhoalpha + rhobeta
      else
         call rho_aop%Gets_dataC_AO(rho_aux(:,:,1))
         call sprepack_ctr('L', M, RMM, rho_aux(MTB+1:MTB+M,MTB+1:MTB+M,1))
      end if

!DFTB: temporary DFTB don't store restarts
      ! Stores the density matrix each 500 steps as a restart.

      call rho_aop%Gets_dataC_ON(rho(:,:,1))
      if (OPEN) call rho_bop%Gets_dataC_ON(rho(:,:,2))

      if ((writedens) .and. ( (mod(istep, td_rst_freq) == 0) .or. &
         (istep.eq.ntdstep) )) then

         if (OPEN) then
            restart_filename='td_a.restart'
            if (istep.eq.ntdstep) restart_filename='td_last_a.restart'
            if (propagator.eq.2) then
               call write_td_restart_magnus(rho(:,:,1), F1a(:,:,1), F1b(:,:,1),&
                                            M, restart_filename)
            else
               call write_td_restart_verlet(rho(:,:,1), M, restart_filename)
            endif
            restart_filename='td_b.restart'
            if (istep.eq.ntdstep) restart_filename='td_last_b.restart'
            if (propagator.eq.2) then
               call write_td_restart_magnus(rho(:,:,2), F1a(:,:,2), F1b(:,:,2),&
                                            M, restart_filename)
            else
               call write_td_restart_verlet(rho(:,:,1), M, restart_filename)
            endif
         else
            restart_filename='td.restart'
            if (istep.eq.ntdstep) restart_filename='td_last.restart'
            if (propagator.eq.2) then
               call write_td_restart_magnus(rho(:,:,1), F1a(:,:,1), F1b(:,:,1),&
                                            M, restart_filename)
            else
               call write_td_restart_verlet(rho(:,:,1), M, restart_filename)
            endif
         end if
      endif

      ! Compute the trace of the density matrix for population analysis.
      if (transport_calc) call transport_rho_trace(M, dim3, rho)

      ! Dipole Moment calculation.
      call td_dipole(t, tdstep, Fx, Fy, Fz, istep, propagator, is_lpfrg, 134)

      ! Population analysis.
      if (transport_calc) call transport_population(M, dim3, natom, Nuc, Iz,   &
                                                    rho_aux, overlap, sqsm,    &
                                                    propagator, is_lpfrg,      &
                                                    istep, OPEN)
      ! TD step finalization.

      if (dftb_calc) call dftb_output(M, dim3,rho_aux, overlap, istep, Iz,     &
                                      natom, Nuc, OPEN)

      call g2g_timer_stop('TD step')
      call g2g_timer_sum_pause("TD - TD Step")

 999  continue

   ! Finalization.
   call write_energies(E1, E2, En, Ens, 0.0D0, Ex, .false., 0.0D0, 0, nsol)

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

subroutine td_allocate_all(M_in,M, dim3, NBCH, propagator, F1a, F1b, fock,     &
                           rho, rho_aux, rhold, rhonew, rho_0, fock_0,sqsm,    &
                           Xmm, Xmat, Xtrans, Ymat, Ytrans, factorial)
   implicit none
   integer, intent(in) :: M, NBCH, propagator, M_in, dim3
   real*8, allocatable, intent(inout) :: F1a(:,:,:), F1b(:,:,:), fock(:,:,:),  &
                                         sqsm(:,:), Xmm(:,:), Xmat(:,:),       &
                                         Xtrans(:,:), Ymat(:,:), Ytrans(:,:),  &
                                         factorial(:), fock_0(:,:,:)
#ifdef TD_SIMPLE
   complex*8 , allocatable, intent(inout) :: rho(:,:,:), rho_aux(:,:,:),       &
                                             rhold(:,:,:), rhonew(:,:,:),      &
                                             rho_0(:,:,:)
#else
   complex*16, allocatable, intent(inout) :: rho(:,:,:), rho_aux(:,:,:),       &
                                             rhold(:,:,:), rhonew(:,:,:),      &
                                             rho_0(:,:,:)
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


   allocate(sqsm(M,M), Xmm(M_in,M_in),Xtrans(M_in,M_in), Xmat(M_in,M_in),      &
            Ymat(M_in,M_in),Ytrans(M_in,M_in), factorial(NBCH))

   allocate(fock(M_in,M_in,dim3), rho(M_in,M_in,dim3),                         &
            rho_aux(M_in,M_in,dim3), rhold(M_in,M_in,dim3),                 &
            rhonew(M_in,M_in,dim3), fock_0(M,M,dim3), rho_0(M,M,dim3))

   if (propagator.eq.2) then
      if ( allocated(F1a) ) deallocate(F1a)
      if ( allocated(F1b) ) deallocate(F1b)

         allocate (F1a(M_in,M_in,dim3), F1b(M_in,M_in,dim3))
   endif

   return
end subroutine td_allocate_all

subroutine td_deallocate_all(F1a, F1b, fock, rho, rho_aux, rhold, &
                             rhonew, Xmm, Xtrans, Ymat, Ytrans,factorial)
   implicit none
   real*8, allocatable, intent(inout) :: F1a(:,:,:), F1b(:,:,:), fock(:,:,:),  &
                                         Xmm(:,:), Xtrans(:,:), Ymat(:,:), &
                                         Ytrans(:,:), factorial(:)
#ifdef TD_SIMPLE
   complex*8 , allocatable, intent(inout) :: rho(:,:,:), rho_aux(:,:,:),    &
                                             rhold(:,:,:), rhonew(:,:,:)
#else
   complex*16, allocatable, intent(inout) :: rho(:,:,:), rho_aux(:,:,:),    &
                                             rhold(:,:,:), rhonew(:,:,:)
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
                          ntatom, natom, Smat, Nuc, a, c, d, Iz, ncont, NORM,  &
                          M, Md, nshell)
    use faint_cpu77, only: intsol
    use faint_cpu  , only: int1
    use mask_ecp   , only: ECP_fock
    implicit none

    double precision, intent(in) :: pc(ntatom), r(ntatom,3)
    integer         , intent(in) :: M, Md, MM, igpu, nsol, ntatom, &
                                    Nuc(M), Iz(natom),nshell(0:4)
    logical         , intent(in) :: NORM
    integer         , intent(inout) :: natom
    double precision, intent(inout) :: RMM11(MM), E1, En, E1s, Ens

    double precision, allocatable, intent(in)    :: a(:,:), c(:,:), d(:,:)
    integer         , allocatable, intent(in)    :: ncont(:)
    double precision, allocatable, intent(inout) :: RMM(:), Smat(:,:)

    integer :: icount

   E1 = 0.0D0 ; En = 0.0D0
   call g2g_timer_sum_start('TD - 1-e Fock')
   call g2g_timer_sum_start('TD - Nuclear attraction')
   call int1(En, RMM, Smat, Nuc, a, c, d, r, Iz, ncont, NORM, natom, M, Md, nshell,ntatom)

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
   return
end subroutine td_get_time

subroutine td_check_prop(is_lpfrg, propagator, istep, lpfrg_steps, tdrestart,&
                         verbose)
   implicit none
   integer, intent(in)  :: propagator, istep, lpfrg_steps, verbose
   logical, intent(in)  :: tdrestart
   logical, intent(out) :: is_lpfrg

   is_lpfrg = ((propagator.eq.1) .or. (((propagator.eq.2) .and. &
              (istep.lt.lpfrg_steps)) .and. (.not.tdrestart)))
   if ( (is_lpfrg).and.(istep.eq.1) ) then
      if (verbose .gt. 2) write(*,'(A)') '  TD - Starting Verlet Propagation'
   endif
   if ( (.not.(is_lpfrg)).and.(((istep-1).eq.lpfrg_steps)) ) then
      if (verbose .gt. 2) write(*,'(A)') '  TD - Starting Magnus Propagation'
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

   ! ELECTRIC FIELD CASE - Perturbation type: Gaussian (default).
   call field_calc(E1, time)

   ! Add 1e contributions to E1.
   do icount = 1, MM
      E1 = E1 + RMM(icount)*RMM11(icount)
   enddo
   E = E1 + E2 + En + Ex
   if (sol) E = E + Es

   return
end subroutine td_calc_energy

subroutine td_dipole(t, tdstep, Fx, Fy, Fz, istep, propagator, is_lpfrg, uid)
   use fileio, only: write_dipole_td, write_dipole_td_header
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
subroutine td_allocate_cublas(M, sizeof_real, sizeof_complex, devPtrX, devPtrY,&
                              devPtrXc, Xmat, Ymat, rho_aux)
   implicit none
   external CUBLAS_INIT, CUBLAS_SHUTDOWN, CUBLAS_ALLOC, CUBLAS_SET_MATRIX
   integer  , intent(in)    :: M, sizeof_real, sizeof_complex
   integer*8, intent(inout) :: devPtrX, devPtrY, devPtrXc
   real*8   , intent(inout) :: Ymat(M,M), Xmat(M,M)
#ifdef TD_SIMPLE
   complex*8 , intent(inout) :: rho_aux(M,M)
#else
   complex*16, intent(inout) :: rho_aux(M,M)
#endif

   integer :: stat, icount, jcount
   stat = 0

   call CUBLAS_INIT()


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

subroutine td_bc_fock_cu(M_in,M, MM, RMM5, fock_op, devPtrX, natom, nshell,    &
                         ncont, istep)
   use dftb_data,        only:dftb_calc, MTB
   use dftb_subs,        only:chimeraDFTB_evol
   use typedef_operator, only:operator
   implicit none
   type(operator), intent(inout) :: fock_op
   integer  , intent(in)      :: M, MM, M_in
   integer*8, intent(in)      :: devPtrX
   real*8   , intent(inout)   :: RMM5(MM)
   real*8 :: fock_0(M,M), fock(M_in,M_in)
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
   call fock_op%Sets_data_AO(fock)
   call fock_op%BChange_AOtoON(devPtrX, M_in,'r')
   call fock_op%Gets_data_ON(fock)
   call sprepack('L', M, RMM5, fock(MTB+1:MTB+M,MTB+1:MTB+M))
   call g2g_timer_stop('fock')

   return
end subroutine td_bc_fock_cu

subroutine td_verlet_cu(M, M_in, dim3, OPEN, fock_aop, rhold, rho_aop, rhonew, &
                        istep, Im, dt_lpfrg, transport_calc, natom, Nuc, Iz,   &
                        overlap, sqsm, devPtrY, devPtrXc, fock_bop, rho_bop)
   use cublasmath       , only : basechange_cublas
   use transport_subs   , only : transport_propagate_cu
   use dftb_data        , only : dftb_calc, rhold_AOTB, rhonew_AOTB
   use typedef_operator , only : operator
   implicit none

   type(operator), intent(inout)           :: fock_aop, rho_aop
   type(operator), intent(inout), optional :: fock_bop, rho_bop

   logical, intent(in)        :: OPEN
   integer, intent(in)        :: dim3
   integer   , intent(in)     :: M,M_in, istep, natom, Nuc(M), Iz(natom)
   integer*8 , intent(in)     :: devPtrY, devPtrXc
   real*8    , intent(in)     :: dt_lpfrg
   logical   , intent(in)     :: transport_calc
   real*8    , intent(inout)  :: overlap(:,:), sqsm(M,M)
#ifdef TD_SIMPLE
   complex*8 , intent(inout)  :: Im, rhold(M_in,M_in, dim3),                   &
                                 rhonew(M_in,M_in, dim3)
   complex*8,  allocatable    :: rho(:,:,:), rho_aux(:,:,:)
#else
   complex*16, intent(inout)  :: Im, rhold(M_in,M_in,dim3),                    &
                                 rhonew(M_in,M_in,dim3)
   complex*16, allocatable    :: rho(:,:,:), rho_aux(:,:,:)
#endif
   integer :: icount, jcount

   allocate(rho(M_in, M_in, dim3), rho_aux(M_in,M_in,dim3))

   call rho_aop%Gets_dataC_ON(rho(:,:,1))
   if(OPEN) call rho_bop%Gets_dataC_ON(rho(:,:,2))

   if(istep.eq.1) then
      call g2g_timer_start('cuconmut')
      call fock_aop%Commut_data_c(rho(:,:,1), rhold(:,:,1), M_in)
      if (OPEN) call fock_bop%Commut_data_c(rho(:,:,2), rhold(:,:,2), M_in)
      call g2g_timer_stop('cuconmut')
      rhold = rho + dt_lpfrg*(Im*rhold)
   endif

   if ((transport_calc) .and. (istep.ge.3)) then
      call rho_aop%Gets_dataC_AO(rho_aux(:,:,1))
      if (OPEN) call rho_bop%Gets_dataC_AO(rho_aux(:,:,2))
      call transport_propagate_cu(M, dim3, natom, Nuc, Iz, 1, istep,           &
                                  overlap, sqsm, rho_aux, devPtrY, OPEN)
   endif

   call g2g_timer_start('commutator')
   call fock_aop%Commut_data_c(rho(:,:,1), rhonew(:,:,1), M_in)
   if (OPEN) call fock_bop%Commut_data_c(rho(:,:,2), rhonew(:,:,2), M_in)
   rhonew = rhold - dt_lpfrg*(Im*rhonew)
   call  g2g_timer_stop('commutator')

   !Transport: Add the driving term to the propagation.
   if ((istep.ge.3) .and. (transport_calc)) then
      write(*,*) 'Transport: Adding driving term to the density.'
      rhonew = rhonew - rho_aux
   endif

   ! Density update (rhold-->rho, rho-->rhonew)
      rhold = rho
      rho   = rhonew

   call rho_aop%Sets_dataC_ON(rho(:,:,1))
   if(OPEN) call rho_bop%Sets_dataC_ON(rho(:,:,2))

!DFTB: rhonew and rhold in AO is store for charge calculations of DFTB
   if (dftb_calc) then
      rhold_AOTB(:,:,1)=basechange_cublas(M_in,rhold(:,:,1),devPtrXc,'inv')
      rhonew_AOTB(:,:,1)=basechange_cublas(M_in,rhonew(:,:,1),devPtrXc,'inv')

      if (OPEN) then
         rhold_AOTB(:,:,2)=basechange_cublas(M_in,rhold(:,:,2),devPtrXc,'inv')
         rhonew_AOTB(:,:,2)=basechange_cublas(M_in,rhonew(:,:,2),devPtrXc,'inv')
      end if
   end if


end subroutine td_verlet_cu

subroutine td_magnus_cu(M, dim3, OPEN,fock_aop, F1a, F1b, rho_aop, rhonew,     &
                        devPtrX, devPtrXc, factorial, NBCH, dt_magnus, natom,  &
                        transport_calc, Nuc, Iz, istep, overlap, sqsm,         &
                        devPtrY, time, M_in, nshell,ncont, fock_bop,           &
                        rho_bop)
   use cublasmath    ,   only: cupredictor, cumagnusfac, basechange_cublas
   use transport_subs,   only: transport_propagate_cu
   use dftb_data,        only: dftb_calc,MTB, rhold_AOTB, rhonew_AOTB
   use dftb_subs,        only: chimeraDFTB_evol
   use typedef_operator, only: operator
   implicit none

   type(operator), intent(inout)           :: fock_aop, rho_aop
   type(operator), intent(inout), optional :: fock_bop, rho_bop

   logical  , intent(in)         :: OPEN
   logical  , intent(in)         :: transport_calc
   integer  , intent(in)         :: M, NBCH, istep, natom, Nuc(M), Iz(natom)
   integer  , intent(in)         :: M_in, dim3
   integer*8, intent(in)         :: devPtrX, devPtrXc, devPtrY
   real*8   , intent(in)         :: dt_magnus, factorial(NBCH), time
   integer  , intent(in)         :: nshell (0:4), ncont(M)
   real*8   , intent(inout)      :: F1a(M_in,M_in,dim3), F1b(M_in,M_in,dim3),  &
                                    overlap(:,:), sqsm(M,M)
#ifdef TD_SIMPLE
   complex*8 , intent(inout) :: rhonew(M_in,M_in,dim3)
   complex*8, allocatable    :: rho(:,:,:), rho_aux(:,:,:)
#else
   complex*16, intent(inout) :: rhonew(M_in,M_in,dim3)
   complex*16, allocatable   :: rho(:,:,:), rho_aux(:,:,:)
#endif
   real*8, allocatable       :: fock_aux(:,:,:), fock(:,:,:)

   allocate(rho(M_in,M_in,dim3), rho_aux(M_in,M_in,dim3),                      &
            fock_aux(M_in,M_in,dim3), fock(M_in,M_in,dim3))

   call fock_aop%Gets_data_ON(fock(:,:,1))
   call rho_aop%Gets_dataC_ON(rho(:,:,1))

   if (OPEN) then
      call fock_bop%Gets_data_ON(fock(:,:,2))
      call rho_bop%Gets_dataC_ON(rho(:,:,2))
   end if

   if (transport_calc) then
      call rho_aop%Gets_dataC_AO(rho_aux(:,:,1))
      if(OPEN) call rho_bop%Gets_dataC_AO(rho_aux(:,:,2))
      call transport_propagate_cu(M, dim3,natom, Nuc, Iz, 2, istep, overlap,   &
                                  sqsm, rho_aux, devPtrY, OPEN)
   endif

! DFTB: this if is temporary, it is to conserve the atomic part of dftb in
! predictor.
   if (dftb_calc) then
      call chimeraDFTB_evol(M,fock(MTB+1:MTB+M,MTB+1:MTB+M,1), fock_aux(:,:,1),&
                            natom, nshell,ncont, istep)
      !DFTB: rhold in AO is store for charge calculations of DFTB
      rhold_AOTB(:,:,1)=basechange_cublas(M_in,rho(:,:,1),devPtrXc,'inv')

      if(OPEN) then
        call chimeraDFTB_evol(M,fock(MTB+1:MTB+M,MTB+1:MTB+M,2),               &
                              fock_aux(:,:,2), natom, nshell,ncont, istep)
        rhold_AOTB(:,:,2)=basechange_cublas(M_in,rho(:,:,2),devPtrXc,'inv')
      end if

      fock=fock_aux
   end if

   call g2g_timer_start('cupredictor')
   call cupredictor(F1a, F1b, fock, rho, devPtrX, factorial, devPtrXc, &
                    dt_magnus, time, M_in, MTB, dim3)
   call g2g_timer_stop('cupredictor')
   call g2g_timer_start('cumagnus')
   call cumagnusfac(fock(:,:,1), rho(:,:,1), rhonew(:,:,1), M_in, NBCH,        &
                    dt_magnus, factorial)
   if(OPEN) call cumagnusfac(fock(:,:,2), rho(:,:,2), rhonew(:,:,2), M_in,     &
                             NBCH, dt_magnus, factorial)
   call g2g_timer_stop('cumagnus')

   !Transport: Add the driving term to the propagation.
   if (transport_calc) then
      write(*,*) 'Transport: Adding driving term to the density.'
      rhonew = rhonew - rho_aux
   endif

!DFTB: rhonew in AO is store for charge calculations of DFTB
   if (dftb_calc) then
      rhonew_AOTB(:,:,1)=basechange_cublas(M_in,rhonew(:,:,1),devPtrXc,'inv')
      if (OPEN) rhonew_AOTB(:,:,2)=basechange_cublas(M_in,rhonew(:,:,2),       &
                                                     devPtrXc,'inv')
   end if

   ! Density update and Fock storage.
   F1a = F1b
   F1b = fock
   rho = rhonew

   call fock_aop%Sets_data_ON(fock(:,:,1))
   call rho_aop%Sets_dataC_ON(rho(:,:,1))

   if (OPEN) then
      call fock_bop%Sets_data_ON(fock(:,:,2))
      call rho_bop%Sets_dataC_ON(rho(:,:,2))
   end if

   return
end subroutine td_magnus_cu

#else

subroutine td_bc_fock(M_in, M, MM, RMM5, fock_op, Xmm, natom, nshell,ncont,    &
                      istep)

   use dftb_data,        only:dftb_calc,MTB
   use dftb_subs,        only:chimeraDFTB_evol
   use typedef_operator, only:operator
   implicit none
   type(operator), intent(inout) :: fock_op
   integer, intent(in)    :: M, MM, M_in
   real*8 , intent(inout) :: RMM5(MM), Xmm(M_in,M_in)
   real*8 , allocatable   :: fock(:,:)
   real*8 , allocatable   :: Xtemp(:,:)
   real*8 , allocatable   :: fock_0(:,:)
   integer, intent(in)  :: natom
   integer, intent(in)  :: ncont(M)
   integer, intent(in)  :: istep
   integer, intent(in)  :: nshell (0:4)

   allocate(fock(M_in,M_in),Xtemp(M_in,M_in),fock_0(M,M))

   call g2g_timer_start('fock')
   call spunpack('L', M, RMM5, fock_0)

   if (dftb_calc) then
      call chimeraDFTB_evol(M,fock_0, fock, natom, nshell,ncont, istep)
   else
      fock=fock_0
   end if

   call fock_op%Sets_data_AO(fock)
   call fock_op%BChange_AOtoON(Xmm, M_in,'r')
   call fock_op%Gets_data_ON(fock)
!carlos: why the ON fock most be store?
   call sprepack('L', M, RMM5, fock(MTB+1:MTB+M,MTB+1:MTB+M))
   call g2g_timer_stop('fock')

   return
end subroutine td_bc_fock

subroutine td_verlet(M, M_in, dim3, OPEN, fock_aop, rhold, rho_aop, rhonew,    &
                     istep, Im, dt_lpfrg, transport_calc, natom, Nuc, Iz,      &
                     overlap, sqsm, Xmat, Ymat, Xtrans, fock_bop, rho_bop)
   use transport_subs,   only : transport_propagate
   use dftb_data     ,   only : dftb_calc, rhold_AOTB, rhonew_AOTB
   use mathsubs      ,   only : basechange_gemm
   use typedef_operator, only: operator
   implicit none

   type(operator), intent(inout) :: fock_aop, rho_aop
   type(operator), intent(inout), optional :: fock_bop, rho_bop

   logical   , intent(in)        :: OPEN
   integer   , intent(in)        :: M, M_in, istep, natom, Nuc(M), Iz(natom),  &
                                    dim3
   real*8    , intent(in)        :: dt_lpfrg
   real*8    , intent(in)        :: Xmat(M_in, M_in), Xtrans(M_in, M_in),      &
                                    Ymat(M_in, M_in)
   logical   , intent(in)        :: transport_calc
   real*8    , intent(inout)     :: overlap(:,:), sqsm(M,M)
#ifdef TD_SIMPLE
   complex*8 , intent(inout)     :: Im, rhold(M_in,M_in,dim3),                 &
                                    rhonew(M_in,M_in,dim3)
   complex*8 , allocatable       :: rho(:,:,:), rho_aux(:,:,:)
#else
   complex*16, intent(inout)     :: Im, rhold(M_in,M_in,dim3),                 &
                                    rhonew(M_in,M_in, dim3)
   complex*16, allocatable       :: rho(:,:,:), rho_aux(:,:,:)
#endif
   integer :: icount, jcount
   real*8 :: traza

   allocate(rho(M_in, M_in, dim3), rho_aux(M_in, M_in, dim3))

!carlos: rho is extracted from rho_op
   call rho_aop%Gets_dataC_ON(rho(:,:,1))

   if(OPEN) call rho_bop%Gets_dataC_ON(rho(:,:,2))

   if (istep.eq.1) then
      call g2g_timer_start('conmutc')
      call fock_aop%Commut_data_c(rho(:,:,1),rhold(:,:,1),M_in)
      if(OPEN) call fock_bop%Commut_data_c(rho(:,:,2),rhold(:,:,2),M_in)
      rhold = rho + dt_lpfrg*(Im*rhold)
      call g2g_timer_stop('conmutc')
   endif

   if ((transport_calc) .and. (istep.ge.3))then
      call rho_aop%Gets_dataC_AO(rho_aux(:,:,1))
      if(OPEN) call rho_bop%Gets_dataC_AO(rho_aux(:,:,2))

      call transport_propagate(M, dim3, natom, Nuc, Iz, 1, istep, overlap,     &
                               sqsm, rho_aux, Ymat, OPEN)
   endif

   call g2g_timer_start('commutator')
   call fock_aop%Commut_data_c(rho(:,:,1),rhonew(:,:,1),M_in)
   if (OPEN) call fock_bop%Commut_data_c(rho(:,:,2),rhonew(:,:,2),M_in)
   rhonew = rhold - dt_lpfrg*(Im*rhonew)
   call  g2g_timer_stop('commutator')

   !Transport: Add the driving term to the propagation.
   if ((istep.ge.3) .and. (transport_calc)) then
      write(*,*) 'Transport: Adding driving term to the density.'
      rhonew = rhonew - rho_aux
   endif
   ! Density update (rhold-->rho, rho-->rhonew)
      rhold = rho
      rho   = rhonew

!carlos: rho_new is now store in rho_op

   call rho_aop%Sets_dataC_ON(rho(:,:,1))
   if(OPEN) call rho_bop%Sets_dataC_ON(rho(:,:,2))

!DFTB: rhonew and rhold in AO is store for charge calculations of DFTB
   if (dftb_calc) then

      rhold_AOTB(:,:,1)=basechange_gemm(M_in,rhold(:,:,1),Xtrans)
      rhonew_AOTB(:,:,1)=basechange_gemm(M_in,rhonew(:,:,1),Xtrans)

      if (OPEN) then
        rhold_AOTB(:,:,1)=basechange_gemm(M_in,rhold(:,:,1),Xtrans)
        rhonew_AOTB(:,:,1)=basechange_gemm(M_in,rhonew(:,:,1),Xtrans)
      end if
   end if

   return
end subroutine td_verlet

subroutine td_magnus(M, dim3, OPEN, fock_aop, F1a, F1b, rho_aop, rhonew,       &
                     factorial, NBCH, dt_magnus, natom, transport_calc, Nuc,   &
                     Iz, istep, overlap, sqsm, Xmat, Ymat, Xtrans, time, M_in, &
                     nshell, ncont, fock_bop, rho_bop)
   use transport_subs,   only: transport_propagate
   use dftb_data,        only:dftb_calc,MTB, rhold_AOTB, rhonew_AOTB
   use dftb_subs,        only:chimeraDFTB_evol
   use mathsubs,         only:basechange_gemm
   use typedef_operator, only:operator
   implicit none

   type(operator), intent(inout)           :: fock_aop, rho_aop
   type(operator), intent(inout), optional :: fock_bop, rho_bop

   logical  , intent(in)      :: transport_calc, OPEN
   integer  , intent(in)      :: dim3
   integer  , intent(in)      :: M, NBCH, istep, natom, Nuc(M), Iz(natom), M_in
   real*8   , intent(in)      :: dt_magnus, factorial(NBCH), Xtrans(M,M), time
   real*8   , intent(in)      :: Xmat(M_in, M_in), Ymat(M_in,M_in)
   real*8   , intent(inout)   :: F1a(M_in,M_in, dim3), F1b(M_in,M_in, dim3),   &
                                 overlap(:,:), sqsm(M,M)
   integer, intent(in)        :: nshell (0:4)
   integer, intent(in)        :: ncont(M)
#ifdef TD_SIMPLE
   complex*8 , intent(inout)  :: rhonew(M_in,M_in,dim3)
   complex*8, allocatable     :: rho(:,:,:), rho_aux(:,:,:), rho_aux0(:,:)
#else
   complex*16, intent(inout)  :: rhonew(M_in,M_in, dim3)
   complex*16, allocatable    :: rho(:,:,:),rho_aux(:,:,:), rho_aux0(:,:)
#endif
   real*8, allocatable        :: fock_aux(:,:,:), fock(:,:,:)
   integer :: ii, jj

   allocate(rho(M_in,M_in,dim3), rho_aux(M_in,M_in,dim3),                      &
            fock_aux(M_in,M_in, dim3), fock(M_in, M_in, dim3),rho_aux0(M_in,M_in))

   call fock_aop%Gets_data_ON(fock(:,:,1))
   call rho_aop%Gets_dataC_ON(rho(:,:,1))

   if (OPEN) then
      call fock_bop%Gets_data_ON(fock(:,:,2))
      call rho_bop%Gets_dataC_ON(rho(:,:,2))
   end if

   if (transport_calc) then
      call rho_aop%Gets_dataC_AO(rho_aux(:,:,1))
      if (OPEN) call rho_bop%Gets_dataC_AO(rho_aux(:,:,2))
      call transport_propagate(M, dim3, natom, Nuc, Iz, 2, istep, overlap,     &
                               sqsm, rho_aux(:,:,1), Ymat, OPEN)
   endif

! DFTB: this if is temporary, it is to conserve the atomic part of dftb in
! predictor.
   if (dftb_calc) then
      call chimeraDFTB_evol(M,fock(MTB+1:MTB+M,MTB+1:MTB+M,1), fock_aux(:,:,1),&
                            natom, nshell,ncont, istep)
!DFTB: rhold in AO is store for charge calculations of DFTB
      rho_aux0=rho(:,:,1)
      rho_aux0=basechange_gemm(M_in,rho_aux0, Xtrans)
      rhold_AOTB(:,:,1) =rho_aux0

      if (OPEN) then
          call chimeraDFTB_evol(M,fock(MTB+1:MTB+M,MTB+1:MTB+M,2),             &
                                fock_aux(:,:,2), natom, nshell,ncont, istep)

        rho_aux0=rho(:,:,2)
        rho_aux0=basechange_gemm(M_in,rho_aux0, Xtrans)
        rhold_AOTB(:,:,2)=rho_aux0
      end if
      fock=fock_aux

   end if

   call g2g_timer_start('predictor')
   call predictor(F1a, F1b, fock, rho, factorial, Xmat, Xtrans, dt_magnus, &
                  time, M_in, MTB, dim3)
   call g2g_timer_stop('predictor')
   call g2g_timer_start('magnus')
   call magnus(fock(:,:,1), rho(:,:,1), rhonew(:,:,1), M_in, NBCH, dt_magnus,  &
               factorial)
   if (OPEN) call magnus(fock(:,:,2), rho(:,:,2), rhonew(:,:,2), M_in, NBCH,   &
                         dt_magnus, factorial)
   call g2g_timer_stop('magnus')

   !Transport: Add the driving term to the propagation.
   if (transport_calc) then
      write(*,*) 'Transport: Adding driving term to the density.'
      rhonew = rhonew - rho_aux
   endif
!DFTB: rhonew in AO is store for charge calculations of DFTB
   if (dftb_calc) then
      rho_aux0 = rhonew(:,:,1)
      rho_aux0=basechange_gemm(M_in,rho_aux0, Xtrans)
      rhonew_AOTB(:,:,1)=rho_aux0
      if (OPEN) then
         rho_aux0 = rhonew(:,:,2)
         rho_aux0=basechange_gemm(M_in,rho_aux0, Xtrans)
         rhonew_AOTB(:,:,2)=rho_aux0
      end if
   end if

   ! Density update and Fock storage.
   F1a = F1b
   F1b = fock
   rho = rhonew
   call fock_aop%Sets_data_ON(fock(:,:,1))
   call rho_aop%Sets_dataC_ON(rho(:,:,1))

   if (OPEN) then
      call fock_bop%Sets_data_ON(fock(:,:,2))
      call rho_bop%Sets_dataC_ON(rho(:,:,2))
   end if

   return
end subroutine td_magnus

#endif

end module time_dependent
