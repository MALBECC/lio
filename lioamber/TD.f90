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
module time_dependent

contains

subroutine TD()
   use garcha_mod    , only: M, Md, NBCH, propagator, tdstep, idip, tdrestart, &
                             exists, RMM, NCO, nang, Iz, natomc, r, d, atmin,  &
                             rmax, jatc, nshell, nnps, nnpp, nnpd, igrid2,     &
                             Nuc, predcoef, npas, nsol, pc, X, Smat, MEMO,     &
                             ntdstep, field, exter, epsilon, writedens, a0,    &
                             sol, kkind, kkinds, cool, cools, GRAD, natom,     &
                             sqsm, Fx, Fy, Fz, Nunp, ntatom
   use ECP_mod       , only: ecpmode, term1e, VAAA, VAAB, VBAC
   use mathsubs
   use transport
   use faint_cpu77   , only: int3lu
   use fileio        , only: read_td_restart_verlet , read_td_restart_magnus , &
                             write_td_restart_verlet, write_td_restart_magnus
#ifdef CUBLAS
   use cublasmath
#endif

   implicit none
   real*8  :: dipxyz(3), q(natom)
   real*8  :: dipole_norm, Qc2, zij, ti, tj, alf, rexp, E, En, E1, E2, E1s,&
              Es, Ens, Ex, ff, t0, t, dt_magnus, dt_lpfrg, tiempo1000, &
              fxx, fyy, fzz, g
   integer :: MM, MMd, M2, M5, M13, M15, M11, unit1, unit2, LWORK,        &
              pert_steps, lpfrg_steps, chkpntF1a, chkpntF1b, igpu,     &
              info, istep, i, j, k, n, ii, jj, kk, jcount, icount
   logical :: ematalloct, dovv, is_lpfrg
   character(len=20) :: restart_filename

   real*8 , allocatable, dimension(:)   :: factorial, Dvec, WORK
   real*8 , allocatable, dimension(:,:) :: Xnano2, Xmm, Xtrans, Ytrans, fock,  &
                                           F1a, F1b, overlap, elmu,&
                                           Ymat, Vmat
! Precision options.
#ifdef TD_SIMPLE
   complex*8  :: Im,Ix
   complex*8 , allocatable, dimension(:,:) :: rho, rho_aux, rhonew, rhold, Xnano
#else
   complex*16 :: Im,Ix
   complex*16, allocatable, dimension(:,:) :: rho, rho_aux, rhonew, rhold, Xnano
#endif

! CUBLAS options.
#ifdef CUBLAS
   integer   :: sizeof_real, sizeof_complex, stat
   integer*8 :: devPtrX, devPtrY, devPtrXc
   external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_SHUTDOWN, CUBLAS_ALLOC,    &
            CUBLAS_GET_MATRIX, CUBLAS_FREE
   parameter(sizeof_real    = 8)
#ifdef TD_SIMPLE
   parameter(sizeof_complex = 8)
#else
   parameter(sizeof_complex = 16)
#endif
#endif

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

   ! Checks and performs allocations.
   call td_allocate_all(M, NBCH, propagator, F1a, F1b, fock, rho, rho_aux,     &
                        rhofirst, rhold, rhonew, sqsm, Vmat, Xmm, Xnano,       &
                        Xnano2, Xtrans, Ymat, Ytrans, Dvec, factorial)
#ifdef CUBLAS
   call td_allocate_cublas(M, sizeof_real, devPtrX, devPtrY)
#endif

   ! Initialises propagator-related parameters and other variables.
   call td_initialise(propagator, tdstep, NBCH, dt_lpfrg, dt_magnus, factorial,&
                      pert_steps, lpfrg_steps, Im, chkpntF1a, chkpntF1b, NCO,  &
                      Nunp, natom, Iz, Qc2)

   ! TD restart reading.
   if (tdrestart) then
      restart_filename = 'td_in.restart'
      if (propagator.eq.2) then
         call read_td_restart_magnus(rho, F1a, F1b, M, restart_filename)
      else
         call read_td_restart_verlet(rho, M, restart_filename)
      endif
      do j = 1, M
      do k = j, M
         if(j.eq.k) then
            RMM(k+(M2-j)*(j-1)/2) = real(rho(j,k))
         else
            RMM(k+(M2-j)*(j-1)/2) = (real(rho(j,k)))*2.0D0
         endif
      enddo
      enddo
   else
      ! Read the density matrix stored in RMM(1,2,3,...,MM) into rho matrix.
      call spunpack_rtc('L', M, RMM, rho)
   endif

   ! Proper TD calculation start.
   write(*,*) 'Starting TD calculation...'

   ! Create integration grid for XC, assigning points to groups (spheres/cubes)
   ! and significant functions to groups, also calculating point weights.
   call td_integration_setup(igrid2, igpu)
   call td_integral_1e(E1, En, E1s, Ens, MM, igpu, nsol, RMM, RMM(M11), r, pc, &
                       ntatom, ecpmode, VAAA, VAAB, VBAC, term1e)

   ! Comment needed here.
   if( transport_calc ) then
      call transport_init(M, natom, Nuc, ngroup, group, mapmat, GammaMagnus,   &
                          GammaVerlet, RMM(M5), overlap, ntdstep)
      call transport_generate_rho(M, rhofirst, rho, generate_rho0)
   endif

   ! Diagonalizes Smat, stores transformation matrix Xmm and calculates the
   ! transposed matrices Xtrans and Ytrans
   call td_overlap_diag(M, Smat, RMM(M13), X, Xtrans, Ymat, Ytrans, Xmm)

   ! Here rho_aux is used as an auxiliar matrix for CUBLAS. Then we need to
   ! restore its original value. Rho is then transformed to the orthonormal
   ! basis with either matmul intrinsic or CUBLAS.
#ifdef CUBLAS
   call td_bc_rho_cu(M, sizeof_real, sizeof_complex, devPtrX, devPtrY, &
                     devPtrXc, X, Ymat, rho, rho_aux)
#else
   call td_bc_rho(M, rho, Ymat, Ytrans)
#endif

   ! Precalculate three-index (two in MO basis, one in density basis) matrix
   ! used in density fitting /Coulomb F element calculation here (t_i in Dunlap)
   call td_coulomb_precalc(igpu, MEMO)

   call g2g_timer_stop('td-inicio')
   ! End of TD initialization.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% TD EVOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   write(*,*) 'TD- Starting propagation.'
   do 999 istep = 1, ntdstep
      ! Checks if step is a leapfrog step.
      is_lpfrg = ((propagator.eq.1) .or. (((propagator.eq.2) .and. &
                 (istep.lt.lpfrg_steps)) .and. (.not.tdrestart)))
      if ( (is_lpfrg).and.(istep.eq.1) ) write(*,*) 'TD - Verlet Propagation'


      call g2g_timer_start('TD step')
      call td_get_time(t, tdstep, istep, propagator, is_lpfrg)

      if (is_lpfrg) then
         call int3lu(E2)
         call g2g_solve_groups(0,Ex,0)
      endif

      ! ELECTRIC FIELD CASE - Perturbation type: Gaussian.
      if ((.not.transport_calc).and.(field)) then
         call td_calc_perturbation(istep, pert_steps, dipxyz, fx, fy, fz,   &
                                   exter, epsilon, a0, Qc2, E1, field, fxx, &
                                   fyy, fzz,g)
      endif

      ! Add 1e contributions to E1.
      do k = 1, MM
         E1 = E1 + RMM(k)*RMM(M11+k-1)
      enddo
      E = E1 + E2 + En + Ex
      if (sol) E = E + Es
      write(*,*) 'TD - Step: ', istep, " Energy : ", E

      ! Here we obtain the Fock matrix in the molecular orbital (MO) basis,
      ! where U matrix with eigenvectors of S, and s is vector with eigenvalues.
      if (is_lpfrg) then
         call g2g_timer_start('fock')
         call spunpack('L', M, RMM(M5), fock)
#ifdef CUBLAS
         Xnano2 = basechange_cublas(M, fock, devPtrX, 'dir')
         fock   = Xnano2
#else
         Xnano2 = matmul(Xtrans, fock)
         fock   = matmul(Xnano2, Xmm)
#endif
         call sprepack('L',M,RMM(M5),fock)
         call g2g_timer_stop('fock')
      endif

      ! Stores F1a and F1b for Magnus propagation.
      if ((propagator.eq.2) .and. (.not.tdrestart)) then
         if (istep.eq.chkpntF1a) F1a = fock
         if (istep.eq.chkpntF1b) F1b = fock
      endif

      if (is_lpfrg) then
         ! Verlet Propagation
         ! In the first step of the propagation we extrapolate rho back in time
         ! using Verlet algorithm to calculate rhold.
         ! Transport: propagation.
#ifdef CUBLAS
         if(istep.eq.1) then
            call g2g_timer_start('cuconmut')
            rhold = commutator_cublas(fock,rho)
            call g2g_timer_stop('cuconmut')
            rhold = rho + dt_lpfrg*(Im*rhold)
         endif

         if ((transport_calc) .and. (istep.ge.3))then
            call transport_propagate_cu(M, natom, Nuc, Iz, ngroup, group,      &
            pop_drive, 1, save_charge_freq, istep, GammaVerlet, overlap, sqsm, &
            rho_aux, rhofirst, devPtrY)
         endif

         call g2g_timer_start('commutator')
         rhonew = commutator_cublas(fock, rho)
         call  g2g_timer_stop('commutator')
#else
         if(istep.eq.1) then
            call g2g_timer_start('conmutc')
            rhold = commutator(fock,rho)
            call g2g_timer_stop('conmutc')
            rhold = rho + dt_lpfrg*(Im*rhold)
         endif
         if ((transport_calc) .and. (istep.ge.3))then
            call transport_propagate(M, natom, Nuc, Iz, ngroup, group,         &
            pop_drive, 1, save_charge_freq, istep, GammaVerlet, overlap, sqsm, &
            rho_aux, rhofirst)
         endif

         call g2g_timer_start('commutator')
         rhonew = commutator(fock, rho)
         call  g2g_timer_stop('commutator')
#endif
         rhonew = rhold - dt_lpfrg*(Im*rhonew)

         !Transport: Add the driving term to the propagation.
         if ((istep.ge.3) .and. (transport_calc)) then
            write(*,*) 'Transport: Adding driving term to the density.'
            rhonew = rhonew - rho_aux
         endif

         ! Density update (rhold-->rho, rho-->rhonew)
         do i = 1, M
         do j = 1, M
            rhold(i,j) = rho(i,j)
            rho(i,j)   = rhonew(i,j)
         enddo
         enddo
      ! End of Verlet propagation.
      else
         ! Magnus Propagation.
         write(*,*) 'TD - Magnus propagation.'

         ! Transport propagation.
         if (transport_calc) then
#ifdef CUBLAS
            call transport_propagate_cu(M, natom, Nuc, Iz, ngroup, group,      &
            pop_drive, 2, save_charge_freq, istep, GammaMagnus, overlap, sqsm, &
            rho_aux, rhofirst, devPtrY)
#else
            call transport_propagate(M, natom, Nuc, Iz, ngroup, group,         &
            pop_drive, 2, save_charge_freq, istep, GammaMagnus, overlap, sqsm, &
            rho_aux, rhofirst)
#endif
         endif

#ifdef CUBLAS
         call g2g_timer_start('cupredictor')
         call cupredictor(F1a, F1b, fock, rho, devPtrX, factorial, fxx, fyy,   &
                          fzz, g, devPtrXc)
         call g2g_timer_stop('cupredictor')
         call g2g_timer_start('cumagnus')
         call cumagnusfac(fock, rho, rhonew, M, NBCH, dt_magnus, factorial)
         call g2g_timer_stop('cumagnus')
#else
         call g2g_timer_start('predictor')
         call predictor(F1a, F1b, fock, rho, factorial, fxx, fyy, fzz, g)
         call g2g_timer_stop('predictor')
         call g2g_timer_start('magnus')
         call magnus(fock, rho, rhonew, M, NBCH, dt_magnus, factorial)
         call g2g_timer_stop('magnus')
#endif
         !Transport: Add the driving term to the propagation.
         if (transport_calc) then
            write(*,*) 'Transport: Adding driving term to the density.'
            rhonew = rhonew - rho_aux
         endif

         ! Density update and Fock storage.
         F1a = F1b
         F1b = fock
         rho = rhonew
      ! End of Magnus Propagation.
      endif

      ! Here we transform the density to the atomic orbital basis and take the
      ! real part of it. The imaginary part of the density can be descarted
      ! since for a basis set of purely real functions the fock matrix is real
      ! and symetric and depends only on the real part of the complex density
      ! matrix. (This will not be true in the case of hybrid functionals)
#ifdef CUBLAS
      call g2g_timer_start('complex_rho_on_to_ao-cu')
      rho_aux = basechange_cublas(M, rho, devPtrXc, 'inv')
      call g2g_timer_stop('complex_rho_on_to_ao-cu')
#else
      call g2g_timer_start('complex_rho_on_to_ao')
      rho_aux = matmul(x(1:M,1:M), rho)
      rho_aux = matmul(rho_aux, Xtrans)
      call g2g_timer_stop('complex_rho_on_to_ao')
#endif
      ! The real part of the density matrix in the atomic orbital basis is
      ! copied in RMM(1,2,3,...,MM) to compute the corresponding Fock matrix.
      do j = 1, M
      do k = j, M
         if (j.eq.k) then
            RMM(k+(M2-j)*(j-1)/2) = real(rho_aux(j,k))
         else
            RMM(k+(M2-j)*(j-1)/2) = (real(rho_aux(j,k)))*2
         endif
      enddo
      enddo

      ! Stores the density matrix each 500 steps as a restart.
      if ((writedens) .and. ((mod(istep,50) == 0) .or. (istep.eq.ntdstep))) then
         restart_filename='td.restart'
         if (istep.eq.ntdstep) restart_filename='td_last.restart'
         if (propagator.eq.2) then
            call write_td_restart_magnus(rho, F1a, F1b, M, restart_filename)
         else
            call write_td_restart_verlet(rho, M, restart_filename)
         endif
      endif

      ! Compute the trace of the density matrix for population analysis.
      if (transport_calc) then
         traza = dcmplx(0.0D0,0.0D0)
         do i = 1, M
            traza = traza + rho(i,i)
         enddo
         write(*,*) 'Trace = ', real(traza)
      endif

      ! Dipole Moment calculation.
      if(istep.eq.1) then
         call write_dipole_td_header(tdstep, Fx, Fy, Fz, 134)
      endif
      if ((propagator.gt.1).and.(is_lpfrg)) then
         if (mod ((istep-1),10) == 0) then
            call g2g_timer_start('DIPOLE_TD')
            call dip(dipxyz)
            call g2g_timer_stop('DIPOLE_TD')
            call write_dipole_td(dipxyz, t, 134)
         endif
      else
         call g2g_timer_start('DIPOLE_TD')
         call dip(dipxyz)
         call g2g_timer_stop('DIPOLE_TD')
         call write_dipole_td(dipxyz, t, 134)
      endif

      ! Population analysis.
      if (transport_calc) then
         if ( ((propagator.gt.1) .and. (is_lpfrg) .and.            &
            (mod((istep-1), save_charge_freq*10) == 0)) .or.                   &
            (mod((istep-1), save_charge_freq) == 0) ) then
            call drive_population(M, natom, Nuc, Iz, pop_drive, ngroup,     &
                                  rho_aux, overlap, group, sqsm)
         endif
      endif

      ! TD step finalization.
      call g2g_timer_stop('TD step')
 999  continue

   ! Finalization.
   call WriteEnergies(E1, E2, En, Ens, 0, Ex, .false., 0)

   deallocate(Xnano, Xnano2, fock, rhonew, rhold, rho, Xmm, Xtrans, Ymat, &
              Ytrans, rho_aux, factorial)
   if (propagator.eq.2) deallocate (F1a, F1b)
#ifdef CUBLAS
   call td_finalise_cublas(devPtrX, devPtrY, devPtrXc)
#endif

   close(134)
   call g2g_timer_stop('TD')

   return
end subroutine TD
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Subroutines used in initialisation.

subroutine td_allocate_all(M, NBCH, propagator, F1a, F1b, fock, rho, rho_aux,  &
                           rhofirst, rhold, rhonew, sqsm, Vmat, Xmm, Xnano,    &
                           Xnano2, Xtrans, Ymat, Ytrans, Dvec, factorial)
   implicit none
   integer, intent(in) :: M, NBCH, propagator
   real*8, allocatable, intent(inout) :: F1a(:,:), F1b(:,:), fock(:,:),        &
                                         sqsm(:,:), Vmat(:,:), Xmm(:,:),       &
                                         Xnano2(:,:), Xtrans(:,:), Ymat(:,:),  &
                                         Ytrans(:,:), Dvec(:), factorial(:)
#ifdef TD_SIMPLE
   complex*8 , allocatable, intent(inout) :: rho(:,:), rho_aux(:,:),    &
                                             rhofirst(:,:), rhold(:,:), &
                                             rhonew(:,:), Xnano(:,:)
#else
   complex*16, allocatable, intent(inout) :: rho(:,:), rho_aux(:,:),    &
                                             rhofirst(:,:), rhold(:,:), &
                                             rhonew(:,:), Xnano(:,:)
#endif

   if ( allocated(fock)     ) deallocate(fock)
   if ( allocated(rho)      ) deallocate(rho)
   if ( allocated(rho_aux)  ) deallocate(rho_aux)
   if ( allocated(rhofirst) ) deallocate(rhofirst)
   if ( allocated(rhold)    ) deallocate(rhold)
   if ( allocated(rhonew)   ) deallocate(rhonew)
   if ( allocated(sqsm)     ) deallocate(sqsm)
   if ( allocated(Vmat)     ) deallocate(Vmat)
   if ( allocated(Xmm)      ) deallocate(Xmm)
   if ( allocated(Xnano)    ) deallocate(Xnano)
   if ( allocated(Xnano2)   ) deallocate(Xnano2)
   if ( allocated(Xtrans)   ) deallocate(Xtrans)
   if ( allocated(Ytrans)   ) deallocate(Ytrans)
   if ( allocated(Ymat)     ) deallocate(Ymat)
   if ( allocated(Dvec)     ) deallocate(Dvec)
   if ( allocated(factorial)) deallocate(factorial)
   allocate(fock(M,M)  , rho(M,M) , rho_aux(M,M) , rhofirst(M,M), rhold(M,M) , &
            rhonew(M,M), sqsm(M,M), Xmm(M,M)     , Xnano(M,M)   , Xnano2(M,M), &
            Xtrans(M,M), Ymat(M,M), Ytrans(M,M)  , Vmat(M,M)    , Dvec(M)    , &
            factorial(NBCH))
   if (propagator.eq.2) then
      if ( allocated(F1a) ) deallocate(F1a)
      if ( allocated(F1b) ) deallocate(F1b)
      allocate (F1a(M,M), F1b(M,M))
   endif

   return
end subroutine td_allocate_all

subroutine td_initialise(propagator, tdstep, NBCH, dt_lpfrg, dt_magnus,        &
                         factorial, pert_steps, lpfrg_steps, Im, chkpntF1a,    &
                         chkpntF1b, NCO, Nunp, natom, Iz, Qc2)
   implicit none
   integer, intent(in)  :: propagator, NBCH, NCO, Nunp, natom, Iz(natom)
   real*8 , intent(in)  :: tdstep
   integer, intent(out) :: pert_steps, lpfrg_steps, chkpntF1a, chkpntF1b
   real*8 , intent(out) :: dt_lpfrg, dt_magnus, factorial(NBCH), Qc2
#ifdef TD_SIMPLE
   complex*8 , intent(out) :: Im
#else
   complex*16, intent(out) :: Im
#endif
   integer :: icount, Nel
   real*8  :: Qc


   ! Common initializations.
   pert_steps  = 100   ; chkpntF1a = 185
   lpfrg_steps = 200   ; chkpntF1b = 195
   Im   = (0.0D0,2.0D0);

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
   Qc2 = 0.0D0

   do icount = 1, natom
      Qc = Qc + Iz(icount)
   enddo
   Qc2 = (Qc - Nel)**2

   return
end subroutine td_initialise

subroutine td_integration_setup(igrid2, igpu)
   implicit none
   integer, intent(in)  :: igrid2
   integer, intent(out) :: igpu

   call g2g_timer_sum_start('Exchange-correlation grid setup')
   call g2g_reload_atom_positions(igrid2)
   call g2g_timer_sum_stop('Exchange-correlation grid setup')

   call aint_query_gpu_level(igpu)
   if (igpu.gt.1) call aint_new_step()

   return
end subroutine td_integration_setup

subroutine td_integral_1e(E1, En, E1s, Ens, MM, igpu, nsol, RMM, RMM11, r, pc, &
                       ntatom, ecpmode, VAAA, VAAB, VBAC, term1e)
   use faint_cpu77, only: int1, intsol

   implicit none
   integer, intent(in)    :: MM, igpu, nsol, ntatom
   logical, intent(in)    :: ecpmode
   real*8 , intent(in)    :: VAAA(MM), VAAB(MM), VBAC(MM), r(ntatom), pc(ntatom)
   real*8 , intent(inout) :: RMM(MM), RMM11(MM), E1, En, E1s, Ens, term1e(MM)
   integer :: icount

   E1 = 0.0D0 ; En = 0.0D0
   call g2g_timer_sum_start('1-e Fock')
   call g2g_timer_sum_start('Nuclear attraction')
   call int1(En)

   ! 1e terms - Pseudopotential terms.
   if (ecpmode) then
      write(*,*) "Adding Pseudopotential terms AAA, AAB, BAC to 1e integrals."
      do icount = 1, MM
         ! Copies 1e terms then adds them the ECP AAA.
         term1e(icount) = RMM11(icount)
         RMM11(icount)  = RMM11(icount) + VAAA(icount) + VAAB(icount) + &
                          VBAC(icount)
      enddo
   end if
   call g2g_timer_sum_stop('Nuclear attraction')

   ! 1e terms - QMMM terms.
   if ((nsol.gt.0) .or. (igpu.ge.4)) then
      call g2g_timer_sum_start('QM/MM')
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
         call g2g_timer_sum_stop('QM/MM')
   endif

   E1=0.D0
   do icount = 1, MM
      E1 = E1 + RMM(icount) * RMM11(icount)
   enddo
   call g2g_timer_sum_stop('1-e Fock')
   return
end subroutine td_integral_1e

subroutine td_overlap_diag(M, Smat, eigenvalues, Xmat, Xtrans, Ymat, Ytrans, Xmm)
   implicit none
   integer, intent(in)    :: M
   real*8 , intent(in)    :: Smat(M,M)
   real*8 , intent(out)   :: Ymat(M,M), Xtrans(M,M), Ytrans(M,M), Xmm(M,M),    &
                             eigenvalues(M)
   real*8 , intent(inout) :: Xmat(M,M)

   integer :: icount, jcount, LWORK, info
   real*8, allocatable :: WORK(:)

   ! Diagonalization of S matrix, both with ESSL or LAPACK.
   ! The S matrix is stored in RMM(M13, M13+1, ..., M13+MM).
   do icount = 1, M
   do jcount = 1, M
      Xmat(icount, jcount) = Smat(icount, jcount)
   enddo
   enddo

   if (allocated(WORK)) deallocate(WORK)
   allocate(WORK(1))
   call dsyev('V', 'L', M, Xmat, M, eigenvalues, WORK, -1, info)

   LWORK = int(WORK(1))
   deallocate(WORK)
   allocate(WORK(LWORK))
   call dsyev('V', 'L', M, Xmat, M, eigenvalues, WORK, LWORK, info)

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
            Xmat(icount, jcount) = 0.0D0
            Ymat(icount, jcount) = 0.0D0
         enddo
      else
         do icount = 1, M
            Ymat(icount,jcount) = Xmat(icount,jcount) *sqrt(eigenvalues(jcount))
            Xmat(icount,jcount) = Xmat(icount,jcount) /sqrt(eigenvalues(jcount))
         enddo
      endif
   enddo

   ! Stores transformation matrix Xmm and transposed matrices.
   do icount = 1, M
   do jcount = 1, M
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

   call g2g_timer_start('Coulomb - precalc')
   if (igpu.gt.2) then
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

subroutine td_calc_perturbation(step, pert_steps, dipxyz, fx, fy, fz, exter, &
                                epsilon, a0, Qc2, E1, field, fxx, fyy, fzz, g)
   use faint_cpu77, only: intfld
   implicit none
   integer, intent(in)    :: step, pert_steps
   logical, intent(inout) :: exter, field
   real*8 , intent(in)    :: epsilon, a0, Qc2
   real*8 , intent(inout) :: dipxyz(3), E1, fx, fy, fz, fxx, fyy, fzz, g
   real*8 :: factor

   E1  = 0.0D0
   fxx = 0.0D0
   fyy = 0.0D0
   fzz = 0.0D0
   if (step.lt.pert_steps) then
      call dip(dipxyz)
      if (exter) then
         g      = 1.0D0
         factor = 2.54D0
         fxx    = fx * exp(-0.2D0*(real(step-50))**2)
         fyy    = fy * exp(-0.2D0*(real(step-50))**2)
         fzz    = fz * exp(-0.2D0*(real(step-50))**2)
         write(*,*) "TD - External field x,y,z in a.u.:"
         write(*,*) fxx, fyy, fzz
      else
         g      = 2.0D0*(epsilon - 1.0D0) / ((2.0D0*epsilon + 1.0D0)*a0**3)
         factor = (2.54D0*2.00D0)
         Fx     = dipxyz(1) / 2.54D0
         Fy     = dipxyz(2) / 2.54D0
         Fz     = dipxyz(3) / 2.54D0
      endif
      call intfld(g, Fxx, Fyy, Fzz)
      E1=-1.00D0 * g * (Fx*dipxyz(1) + Fy*dipxyz(2) + Fz*dipxyz(3)) / factor - &
          0.50D0 * (1.0D0 - 1.0D0/epsilon) * Qc2/a0
   else
      field = .false.
   endif

end subroutine td_calc_perturbation

! CUBLAs-dependent subroutines.
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


subroutine td_finalise_cublas(devPtrX, devPtrY, devPtrXc)
   implicit none
   external CUBLAS_SHUTDOWN, CUBLAS_FREE
   integer*8 :: devPtrX, devPtrXc, devPtrY
   call CUBLAS_FREE(devPtrX)
   call CUBLAS_FREE(devPtrXc)
   call CUBLAS_FREE(devPtrY)
   call CUBLAS_SHUTDOWN()
end subroutine td_finalise_cublas
#endif

end module time_dependent
