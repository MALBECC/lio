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
                             sqsm, Fx, Fy, Fz, Nunp
   use ECP_mod       , only: ecpmode, term1e, VAAA, VAAB, VBAC
   use mathsubs
   use transport
   use faint_cpu77   , only: int1, intsol, int3mem, intfld, int3lu
   use fileio        , only: read_td_restart_verlet , read_td_restart_magnus , &
                             write_td_restart_verlet, write_td_restart_magnus
#ifdef CUBLAS
   use cublasmath
#endif

   implicit none
   real*8  :: dipxyz(3), q(natom)
   real*8  :: dipole_norm, Qc2, zij, ti, tj, alf, rexp, E, En, E1, E2, E1s,&
              Es, Ens, Ex, g, factor, fxx, fyy, fzz, ff, t0, t, dt_magnus,     &
              dt_lpfrg, tiempo1000
   integer :: MM, MMd, M2, M5, M13, M15, M11, unit1, unit2, LWORK,        &
              pert_steps, lpfrg_steps, chkpntF1a, chkpntF1b, igpu,     &
              info, istep, i, j, k, n, ii, jj, kk
   logical :: ematalloct, dovv
   character(len=20) :: restart_filename

   real*8 , allocatable, dimension(:)   :: factorial, Dvec, WORK
   real*8 , allocatable, dimension(:,:) :: Xnano2, Xmm, Xtrans, Ytrans, fock,  &
                                           F1a, F1b, overlap, elmu,&
                                           Ymat, Vmat

! Precision options.
#ifdef TD_SIMPLE
   complex*8  :: Im,Ix
   complex*8 , allocatable, dimension(:,:) :: rho, rho1, rhonew, rhold, Xnano
#else
   complex*16 :: Im,Ix
   complex*16, allocatable, dimension(:,:) :: rho, rho1, rhonew, rhold, Xnano
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
   call td_allocate_all(M, NBCH, propagator, F1a, F1b, fock, rho, rho1,        &
                        rhofirst, rhold, rhonew, sqsm, Vmat, Xmm, Xnano,       &
                        Xnano2, Xtrans, Ymat, Ytrans, Dvec, factorial)
#ifdef CUBLAS
   call td_allocate_cublas(M, sizeof_real, devPtrX, devPtrY)
#endif

   ! Initializations.
   pert_steps  = 100   ; chkpntF1a = 185
   lpfrg_steps = 200   ; chkpntF1b = 195
   E           = 0.0D0 ; E1        = 0.0D0
   En          = 0.0D0 ; E2        = 0.0D0
   Im   = (0.0D0,2.0D0)

   ! Initialises propagator-related parameters.
   call init_propagators(propagator, tdstep, NBCH, dt_lpfrg, dt_magnus, factorial)

   ! Gets squared total atomic charge.
   call get_atomic_charge(NCO, Nunp, natom, Iz, Qc2)

   ! Transport Initializations
   if (transport_calc) call transport_init(M, natom, Nuc, ngroup, group, mapmat, &
                                           GammaMagnus, GammaVerlet)

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

   ! Two electron integral with neighbor list.
   do i = 1, natom
      natomc(i) = 0
      do j = 1, natom
         d(i,j) = (r(i,1)-r(j,1))**2 + (r(i,2)-r(j,2))**2 + (r(i,3)-r(j,3))**2
         zij    = atmin(i) + atmin(j)
         ti     = atmin(i) / zij
         tj     = atmin(j) / zij
         alf    = atmin(i) * tj
         rexp   = alf      * d(i,j)
         if (rexp.lt.rmax) then
            natomc(i)         = natomc(i) + 1
            jatc(natomc(i),i) = j
         endif
      enddo
   enddo

   ! Basis index reorganization.
   do ii = nshell(0), 1, -1
      nnps(Nuc(ii)) = ii
   enddo
   do ii = nshell(0) + nshell(1), nshell(0) + 1, -1
      nnpp(Nuc(ii)) = ii
   enddo
   do ii = M, nshell(0) + nshell(1) + 1, -1
      nnpd(Nuc(ii)) = ii
   enddo

   ! Create integration grid for XC, assigning points to groups (spheres/cubes)
   ! and significant functions to groups, also calculating point weights.
   call g2g_timer_sum_start('Exchange-correlation grid setup')
   call g2g_reload_atom_positions(igrid2)
   call g2g_timer_sum_stop('Exchange-correlation grid setup')

   call aint_query_gpu_level(igpu)
   if (igpu.gt.1) call aint_new_step()

   ! 1e terms - HH core, 1 electron matrix elements
   call g2g_timer_sum_start('1-e Fock')
   call g2g_timer_sum_start('Nuclear attraction')
   call int1(En)

   ! 1e terms - Pseudopotential terms.
   if (ecpmode) then
      write(*,*) "Adding Pseudopotential terms AAA and AAB to 1e integrals."
      do k=1,MM
         ! Copies 1e termsm then adds them the ECP AAA.
         term1e(k)    = RMM(M11+k-1)
         RMM(M11+k-1) = RMM(M11+k-1) + VAAA(k) + VAAB(k) + VBAC(k)
      enddo
   end if

   ! 1e terms - QMMM terms.
   call g2g_timer_sum_stop('Nuclear attraction')
   if(nsol.gt.0.or.igpu.ge.4) then
      call g2g_timer_sum_start('QM/MM')
      if (igpu.le.1) then
         call g2g_timer_start('intsol')
         call intsol(E1s,Ens,.true.)
         call g2g_timer_stop('intsol')
      else
         call aint_qmmm_init(nsol,r,pc)
         call g2g_timer_start('aint_qmmm_fock')
         call aint_qmmm_fock(E1s,Ens)
         call g2g_timer_stop('aint_qmmm_fock')
      endif
         call g2g_timer_sum_stop('QM/MM')
   endif

   E1=0.D0
   do k = 1, MM
      E1 = E1 + RMM(k)*RMM(M11+k-1)
   enddo
   call g2g_timer_sum_stop('1-e Fock')

   ! Comment needed here.
   if( transport_calc )then
      if (allocated(overlap)) deallocate(overlap)
      allocate(overlap(M,M))
      call spunpack('L', M, RMM(M5), overlap)
   endif

   ! Diagonalization of S matrix, both with ESSL or LAPACK.
   ! The S matrix is stored in RMM(M13, M13+1, ..., M13+MM).
#ifdef essl
   call DSPEV(1, RMM(M5), RMM(M13), X, M, M, RMM(M15), M2)
#endif
#ifdef pack
   do ii = 1, M
   do jj = 1, M
      X(ii,jj) = Smat(ii,jj)
   enddo
   enddo

   if (allocated(WORK)) deallocate(WORK)
   allocate(WORK(1))
   call dsyev('V', 'L', M, X, M, RMM(M13), WORK, -1, info)

   LWORK = int(WORK(1))
   deallocate(WORK)
   allocate(WORK(LWORK))
   call dsyev('V', 'L', M, X, M, RMM(M13), WORK, LWORK, info)
#endif

   ! Here we obtain the transformation matrices X and Y for converting from the
   ! atomic orbital basis to the molecular orbital basis (truncated during
   ! linear dependency elimination). S is the overlap matrix, s is the diagonal
   ! eigenvalue matrix of S and U is the eigenvector matrix of S:
   ! X = U s^(-1/2)
   ! Matrix X's dimension is M*3M. In the first M*M terms it contains the
   ! transformation matrices and in the other M*2M terms it contains auxiliar
   ! matrices.
   call g2g_timer_start('inicio1')
   do j = 1, M
      if (RMM(M13+j-1).lt.1.0D-06) then
         write(*,*) 'WARNING - TD: Linear dependency detected in S matrix.'
         do i = 1, M
            X(i,j) = 0.0D0
            Ymat(i,j) = 0.0D0
         enddo
      else
         do i = 1, M
            Ymat(i,j) = X(i,j) * sqrt(RMM(M13+j-1))
            X(i,j) = X(i,j) / sqrt(RMM(M13+j-1))
         enddo
      endif
   enddo

   ! Here rho1 is used as an auxiliar matrix. Then we need to restore its
   ! original value.
#ifdef CUBLAS
   do i = 1, M
   do j = 1, M
      rho1(i,j) = cmplx(X(i,j), 0.0D0)
   enddo
   enddo

   call CUBLAS_ALLOC(M*M, sizeof_real, devPtrX)
   call CUBLAS_ALLOC(M*M, sizeof_complex, devPtrXc)
   call CUBLAS_ALLOC(M*M, sizeof_complex, devPtrY)
   if (stat.NE.0) then
      write(*,*) "ERROR - TD: X and/or Y CUBLAS memory allocation failed."
      call CUBLAS_SHUTDOWN()
      stop
   endif

   call CUBLAS_SET_MATRIX(M, M, sizeof_complex, rho1, M, devPtrXc, M)
   call CUBLAS_SET_MATRIX(M, M, sizeof_real, x, M, devPtrX, M)
   do i = 1, M
   do j = 1, M
      rho1(i,j) = cmplx(Ymat(i,j), 0.0D0)
   enddo
   enddo

   call CUBLAS_SET_MATRIX(M, M, sizeof_complex, rho1, M, devPtrY, M)
   if (stat.NE.0) then
      write(*,*) "ERROR - TD : X and/or Y CUBLAS setting failed."
      call CUBLAS_SHUTDOWN
      stop
   endif
   rho1 = 0D0
#endif

   ! The transformation matrix is stored in Xmm, and the transposed matrices are
   ! calculated.
   do i = 1, M
   do j = 1, M
      Xmm(i,j)    = X(i,j)
      Xtrans(j,i) = X(i,j)
      Ytrans(j,i) = Ymat(i,j)
   enddo
   enddo

   ! Steps needed for transport.
   if (transport_calc) call transport_generate_rho(M, rhofirst, rho, generate_rho0)

   ! Rho is transformed to the orthonormal basis with either matmul intrinsic
   ! or CUBLAS:
#ifdef CUBLAS
   call g2g_timer_start('complex_rho_ao_to_on-cu')
   rho1 = basechange_cublas(M, rho, devPtrY, 'dir')
   rho  = rho1
   call g2g_timer_stop('complex_rho_ao_to_on-cu')
#else
   rho = matmul(Ytrans, rho)
   rho = matmul(rho   , Ymat)
#endif

   ! Precalculate three-index (two in MO basis, one in density basis) matrix
   ! used in density fitting /Coulomb F element calculation here (t_i in Dunlap)
   call aint_query_gpu_level(igpu)
   if (igpu.gt.2) then
      call aint_coulomb_init()
      if (igpu.eq.5) MEMO = .false.
   endif
   if (MEMO) then
      call g2g_timer_start('int3mem')
      call g2g_timer_sum_start('Coulomb precalc')
      call int3mem()
      call g2g_timer_stop('int3mem')
      call g2g_timer_sum_stop('Coulomb precalc')
   endif

#ifdef CUBLAS
   if (.not.transport_calc) call CUBLAS_FREE(devPtrY)
#endif

   call g2g_timer_stop('td-inicio')
   ! End of TD initialization.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% TD EVOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   write(*,*) 'TD- Strarting propagation.'
   do 999 istep=1, ntdstep

      call g2g_timer_start('TD step')
      if ((propagator.eq.2).and.(istep.lt.lpfrg_steps).and.(.not.tdrestart))   &
      then
         t = (istep-1) * tdstep * 0.1
      else
         t = 20 * tdstep
         t = t + (istep-200)*tdstep
      endif
      if (propagator.eq.1) then
         t = (istep-1) * tdstep
      endif
      t = t * 0.024190D0
      write(*,*) 'TD - Time (fs)  =', t

      if ((propagator.eq.1).or. (((propagator.eq.2).and.(istep.lt.lpfrg_steps))&
         .and. (.not.tdrestart))) then
         call int3lu(E2)
         call g2g_solve_groups(0,Ex,0)
      endif

      write(*,*) 'TD - Step: ', istep, " Energy : ", E+Ex
      E1 = 0.0D0

      ! ELECTRIC FIELD CASE - Perturbation type: Gaussian.
      fxx = 0.0D0
      fyy = 0.0D0
      fzz = 0.0D0
      if(.not.transport_calc) then
      if(istep.lt.pert_steps) then
         if (field) then
            call dip(dipxyz)
            if (exter) then
               g      = 1.0D0
               factor = 2.54D0
               fxx    = fx * exp(-0.2D0*(real(istep-50))**2)
               fyy    = fy * exp(-0.2D0*(real(istep-50))**2)
               fzz    = fz * exp(-0.2D0*(real(istep-50))**2)
               write(*,*) "TD - External field x,y,z in a.u.: ", fxx, fyy, fzz
            else
               g  = 2.0D0 * (epsilon-1.0D0) / ((2.0D0*epsilon+1.0D0)*a0**3)
               Fx = dipxyz(1) / 2.54D0
               Fy = dipxyz(2) / 2.54D0
               Fz = dipxyz(3) / 2.54D0
               factor = (2.54D0*2.00D0)
            endif
            write(*,*) 'TD - Epsilon =', epsilon
            call intfld(g, Fxx, Fyy, Fzz)
            !TO-DO: shouldn't E1 use Fxx instead of Fx?? (idem y, z)
            E1 = -1.00D0 * g * (Fx*dipxyz(1) + Fy*dipxyz(2) + Fz*dipxyz(3)) &
                  / factor - 0.50D0 * (1.0D0-1.0D0/epsilon) * Qc2 / a0
         endif
      else
         field = .false.
      endif
      endif

      ! E1 includes solvent 1e contributions.
      do k = 1, MM
         E1 = E1 + RMM(k)*RMM(M11+k-1)
      enddo

      ! Here we obtain the Fock matrix in the molecular orbital (MO) basis,
      ! where U matrix with eigenvectors of S, and s is vector with eigenvalues.
      if ((propagator.eq.1).or. (((propagator.eq.2).and.(istep.lt.lpfrg_steps))&
         .and. (.not.tdrestart))) then
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

      ! Fock is stored in molecular orbital basis. Also, stores F1a and F1b
      ! for Magnus propagation.
      if (propagator.eq.2) then
         if (.not.tdrestart) then
            if (istep.eq.chkpntF1a) then
               F1a = fock
            endif
            if (istep.eq.chkpntF1b) then
               F1b = fock
            endif
         endif

         ! Writes a checkpoint restart
         if ((writedens).and.(mod (istep, 500) == 0)) then
            call write_td_restart_magnus(rho, F1a, F1b, M, restart_filename)
         endif
      endif

      E = E1 + E2 + En
      if (sol) E = E + Es

      if ((propagator.eq.1).or. (((propagator.eq.2).and.(istep.lt.lpfrg_steps))&
         .and.(.not.tdrestart))) then
         ! Verlet Propagation
         write(*,*) 'TD - Verlet Propagation'

         ! In the first step of the propagation we extrapolate rho back in time
         ! using Verlet algorithm to calculate rhold.
         if(istep.eq.1) then
            rhold = rho + dt_lpfrg*(Im*rhold)
#ifdef CUBLAS
            call g2g_timer_start('cuconmut')
            rhold = commutator_cublas(fock,rho)
            call g2g_timer_stop('cuconmut')
#else
            call g2g_timer_start('conmutc')
            rhold = commutator(fock,rho)
            rhold = rho + dt_lpfrg*(Im*rhold)
            call g2g_timer_stop('conmutc')
#endif
         endif

         ! Transport: propagation.
         if (transport_calc) then
            call g2g_timer_start('TRANSPORT - b Verlet -')
            if(istep.eq.1) then
               unit1 = 55555
               if ( pop_drive == 1 ) then
                  open( unit = unit1, file = 'DriveMul')
               elseif ( pop_drive == 2 ) then
                  open(unit = unit1, file = 'DriveLowd')
               endif
            endif
            if(istep.ge.3) then
               scratchgamma = GammaVerlet * exp(-0.0001D0*(dble(istep-1000))**2)
               call ELECTROSTAT(rho1, mapmat, overlap, rhofirst, scratchgamma,M)
               if ( (pop_drive == 1) .or. (pop_drive == 2) ) then
                  if (mod(istep-1,save_charge_freq*10) == 0) then
                     q = 0.0D0
                     call drive_population(M, natom, Nuc, Iz, pop_drive, ngroup, &
                                           rho1, overlap, group, sqsm)
                     endif
               endif
            endif

#ifdef CUBLAS
            call g2g_timer_start('complex_rho_ao_to_on-cu')
            rho1 = basechange_cublas(M, rho1, devPtrY, 'dir')
            call g2g_timer_stop('complex_rho_ao_to_on-cu')
#endif
            call g2g_timer_stop('TRANSPORT - b Verlet -')
         endif

         ! Using commutator:
         call g2g_timer_start('commutator')
#ifdef CUBLAS
         rhonew = commutator_cublas(fock,rho)
         rhonew = rhold - dt_lpfrg*(Im*rhonew)
#else
         rhonew = commutator(fock,rho)
         rhonew = rhold - dt_lpfrg*(Im*rhonew)
#endif
         call  g2g_timer_stop('commutator')

         !Transport: Add the driving term to the propagation.
         if ((istep.ge.3) .and. (transport_calc)) then
            write(*,*) 'Transport: Adding driving term to the density.'
            rhonew = rhonew - rho1
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
            call g2g_timer_start('TRANSPORT - b magnus -')
            if(istep.le.1000) then
               scratchgamma = GammaMagnus*exp(-0.0001D0*(dble(istep-1000))**2)
            else
               scratchgamma = GammaMagnus
            endif
            call ELECTROSTAT(rho1, mapmat, overlap, rhofirst, scratchgamma, M)
            if (pop_drive == 1 .or. pop_drive == 2) then
               if (mod(istep-1,save_charge_freq)==0) then
                  call drive_population(M, natom, Nuc, Iz, pop_drive, ngroup, &
                                        rho1, overlap, group, sqsm)
               endif
            endif
#ifdef CUBLAS
            call g2g_timer_start('complex_rho_ao_to_on-cu')
            rho1 = basechange_cublas(M, rho1, devPtrY, 'dir')
            call g2g_timer_stop('complex_rho_ao_to_on-cu')
#endif
            call g2g_timer_stop('TRANSPORT - b magnus -')
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
            rhonew = rhonew - rho1
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
      rho1 = basechange_cublas(M, rho, devPtrXc, 'inv')
      call g2g_timer_stop('complex_rho_on_to_ao-cu')
#else
      call g2g_timer_start('complex_rho_on_to_ao')
      rho1 = matmul(x(1:M,1:M), rho)
      rho1 = matmul(rho1, Xtrans)
      call g2g_timer_stop('complex_rho_on_to_ao')
#endif
      ! The real part of the density matrix in the atomic orbital basis is
      ! copied in RMM(1,2,3,...,MM) to compute the corresponding Fock matrix.
      do j = 1, M
      do k = j, M
         if (j.eq.k) then
            RMM(k+(M2-j)*(j-1)/2) = real(rho1(j,k))
         else
            RMM(k+(M2-j)*(j-1)/2) = (real(rho1(j,k)))*2
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
      if ((propagator.eq.2) .and. (istep.lt.lpfrg_steps) .and.(.not.tdrestart))&
      then
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
         if ( ((propagator.eq.2) .and. (istep.lt.lpfrg_steps) .and. (.not.tdrestart) .and.&
         (mod((istep-1), save_charge_freq*10) == 0)) .or.                                 &
         (mod((istep-1), save_charge_freq) == 0) ) then
            call drive_population(M, natom, Nuc, Iz, pop_drive, ngroup, rho1, overlap,    &
                                  group, sqsm)
         endif
      endif

      ! TD step finalization.
      call g2g_timer_stop('TD step')
      write(*,*)
      if (istep .eq. 1000) then
         call g2g_timer_start('corrida 1000')
         tiempo1000 = t
      elseif (istep .eq. 2000) then
         call g2g_timer_stop('corrida 1000')
         write(*,*) t - tiempo1000
      endif

      ! Transport: stops TD and indicates the generation of RHO0
      if((istep.ge.1).and.(generate_rho0)) then
         write(*,*) "Transport - RHO0 generated."
         exit
      endif

 999  continue

   ! Finalization.
   if (memo) deallocate(kkind, kkinds, cool ,cools)
   if (propagator.eq.2) deallocate (F1a, F1b)

   if (GRAD) then
      write(*,*)
      write(*,600)
      write(*,610)
      write(*,620) E1, E2-Ex, En
      if (sol) then
         write(*,615)
         write(*,625) Es
      endif
      write(*,*)
      write(*,450) E
   else
      E = E - Ex
   endif

   ! Calculation of energy weighted density matrix,
   kk =0
   do j = 1, M
   do i = j, M
      kk = kk+1
      RMM(M15+kk-1) = 0.D0
      if(i.eq.j) then
         ff = 2.D0
      else
         ff = 4.D0
      endif
      do k = 1, NCO
         RMM(M15+kk-1) = RMM(M15+kk-1) - RMM(M13+k-1)*ff*X(i,M2+k)*X(j,M2+k)
      enddo
   enddo
   enddo

#ifdef CUBLAS
   call CUBLAS_FREE(devPtrX)
   call CUBLAS_FREE(devPtrXc)
   call CUBLAS_FREE(devPtrY)
   call CUBLAS_SHUTDOWN()
#endif

   close(134)
   call g2g_timer_stop('TD')
   deallocate(Xnano, Xnano2, fock, rhonew, rhold, rho, Xmm, Xtrans, Ymat, &
              Ytrans, rho1, factorial)
   return

! Formats for Energy Printing
 450  format ('FINAL ENERGY = ',F19.12)
 600  format('  ENERGY CONTRIBUTIONS IN A.U.')
 610  format(2x,'ONE ELECTRON',9x,'COULOMB',11x,'NUCLEAR')
 615  format(2x,'SOLVENT')
 620  format(F14.7,4x,F14.7,4x,F14.7)
 625  format(F14.7)
end subroutine TD
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Checks and performs allocations.
subroutine td_allocate_all(M, NBCH, propagator, F1a, F1b, fock, rho, rho1,     &
                           rhofirst, rhold, rhonew, sqsm, Vmat, Xmm, Xnano,    &
                           Xnano2, Xtrans, Ymat, Ytrans, Dvec, factorial)
   implicit none
   integer, intent(in) :: M, NBCH, propagator
   real*8, allocatable, intent(inout) :: F1a(:,:), F1b(:,:), fock(:,:),        &
                                         sqsm(:,:), Vmat(:,:), Xmm(:,:),       &
                                         Xnano2(:,:), Xtrans(:,:), Ymat(:,:),  &
                                         Ytrans(:,:), Dvec(:), factorial(:)
#ifdef TD_SIMPLE
   complex*8 , allocatable, intent(inout) ::rho(:,:), rho1(:,:), rhofirst(:,:),&
                                            rhold(:,:), rhonew(:,:), Xnano(:,:)
#else
   complex*16, allocatable, intent(inout) ::rho(:,:), rho1(:,:), rhofirst(:,:),&
                                            rhold(:,:), rhonew(:,:), Xnano(:,:)
#endif

   if ( allocated(fock)     ) deallocate(fock)
   if ( allocated(rho)      ) deallocate(rho)
   if ( allocated(rho1)     ) deallocate(rho1)
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
   allocate(fock(M,M)  , rho(M,M) , rho1(M,M)  , rhofirst(M,M), rhold(M,M) ,   &
            rhonew(M,M), sqsm(M,M), Xmm(M,M)   , Xnano(M,M)   , Xnano2(M,M),   &
            Xtrans(M,M), Ymat(M,M), Ytrans(M,M), Vmat(M,M)    , Dvec(M)    ,   &
            factorial(NBCH))
   if (propagator.eq.2) then
      if ( allocated(F1a) ) deallocate(F1a)
      if ( allocated(F1b) ) deallocate(F1b)
      allocate (F1a(M,M), F1b(M,M))
   endif

   return
end subroutine td_allocate_all


subroutine init_propagators(propagator, tdstep, NBCH, dt_lpfrg, dt_magnus, factorial)
   implicit none
   integer, intent(in)  :: propagator, NBCH
   real*8 , intent(in)  :: tdstep
   real*8 , intent(out) :: dt_lpfrg, dt_magnus, factorial(NBCH)

   integer :: icount

   ! Initialises propagator-related parameters.
   select case (propagator)
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
         write(*,*) "ERROR - TD: Wrong value for propagator (init_propagators)."
         stop
   end select

   return
end subroutine init_propagators

subroutine get_atomic_charge(NCO, Nunp, natom, Iz, Qc2)
   implicit none
   integer, intent(in)  :: NCO, Nunp, natom, Iz(natom)
   real*8 , intent(out) :: Qc2
   integer :: Nel, i
   real*8  :: Qc

   Nel = 2*NCO + Nunp
   Qc  = 0.0D0
   Qc2 = 0.0D0

   do i = 1, natom
      Qc = Qc + Iz(i)
   enddo
   Qc2 = (Qc - Nel)**2

   return
end subroutine get_atomic_charge

#ifdef CUBLAS
subroutine td_allocate_cublas(M, sizeof_real, devPtrX, devPtrY)

   implicit none
   integer, intent(in)    :: M, sizeof_real
   integer, intent(inout) :: devPtrX, devPtrY
   integer :: stat
   stat = 0

   write(*,*) 'TD: Using CUBLAS.'
   call CUBLAS_INIT()
   call CUBLAS_ALLOC(M*M, sizeof_real, devPtrX)
   call CUBLAS_ALLOC(M*M, sizeof_real, devPtrY)
   if (stat.NE.0) then
      write(*,*) "ERROR - TD: CUBLAS initialization failed."
      call CUBLAS_SHUTDOWN
      stop
   endif

   return
end subroutine td_allocate_cublas
#endif

end module time_dependent
