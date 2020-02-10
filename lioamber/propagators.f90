!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% MAGNUS.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains the routine for an N-order magnus propagation.            !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module propagators
contains

subroutine magnus(fock, rhoOld, rhoNew, M, N, dt, factorial)
    ! Input:  Fock(t+(deltat/2)), rho(t)
    ! Output: rho6 = rho(t+deltat)
    implicit none
    integer         , intent(in)  :: M, N
    double precision, intent(in)  :: Fock(M,M), dt, factorial(N)
#ifdef TD_SIMPLE
    complex         , intent(in)  :: RhoOld(M,M)
    complex         , intent(out) :: RhoNew(M,M)
    complex, parameter   :: ICMPLX = CMPLX(0.0D0,1.0D0)
    complex, allocatable :: Scratch(:,:), ConmPrev(:,:), ConmNext(:,:), &
                            Omega1(:,:)
#else
    double complex  , intent(in)  :: RhoOld(M,M)
    double complex  , intent(out) :: RhoNew(M,M)
    double complex, parameter   :: ICMPLX = CMPLX(0.0D0,1.0D0)
    double complex, allocatable :: Scratch(:,:), ConmPrev(:,:), ConmNext(:,:), &
                                   Omega1(:,:)
#endif
    integer :: icount

    ! Variable initializations.
    allocate(ConmPrev(M,M), ConmNext(M,M), Omega1(M,M), Scratch(M,M))
    ConmPrev = RhoOld
    RhoNew   = RhoOld

    ! Omega1 (=W) instantiation.
    Omega1 = -1 * ICMPLX * Fock * dt

    ! Density matrix propagation
    do icount = 1, N
        ConmNext = MATMUL(Omega1, ConmPrev)
        Scratch  = MATMUL(ConmPrev, Omega1)
        ConmNext = ConmNext - Scratch
        RhoNew   = RhoNew + factorial(icount) * ConmNext
        ConmPrev = ConmNext
    enddo

    deallocate(ConmPrev, ConmNext, Omega1, Scratch)
    return
end subroutine magnus

subroutine predictor(F1a, F1b, FON, rho2, factorial, Xmat, Xtrans, timestep, &
                     time, M_in, MTB, dim3)
   ! This routine receives: F1a, F1b, rho2
   ! And gives: F5 = F(t+(deltat/2))
   use garcha_mod   , only: NBCH, rhoalpha, rhobeta, OPEN, r, d, natom,      &
                            ntatom, Iz, MEMO, Fmat_vec, Fmat_vec2, Ginv_vec, &
                            Hmat_vec, Gmat_vec, Pmat_vec
   use field_subs   , only: field_calc
   use mathsubs     , only: basechange
   use faint_cpu    , only: int3lu
   use fockbias_subs, only: fockbias_apply
   use basis_data   , only: M

   implicit none
   ! MTB is only used in DFTB, it equals 0 otherwise.
   integer         , intent(in)    :: M_in, dim3, MTB
   double precision, intent(in)    :: Xtrans(M_in,M_in), timestep, time, &
                                      factorial(NBCH), Xmat(M_in,M_in)
   double precision, intent(inout) :: F1a(M_in,M_in,dim3), F1b(M_in,M_in,dim3),&
                                      FON(M_in,M_in,dim3)
#ifdef TD_SIMPLE
   complex, intent(in)  :: rho2(M_in,M_in,dim3)
   complex, allocatable :: rho4(:,:,:), rho2t(:,:,:)
#else
   double complex, intent(in)  :: rho2(M_in,M_in,dim3)
   double complex, allocatable :: rho4(:,:,:), rho2t(:,:,:)
#endif
   integer :: M2
   double precision :: E2, tdstep1, Ex, E1
   double precision, allocatable :: F3(:,:,:), FBA(:,:,:)

   allocate(rho4(M_in,M_in,dim3), rho2t(M_in,M_in,dim3), F3(M_in,M_in,dim3), &
            FBA(M_in,M_in,dim3))
   M2 = 2 * M

   ! Initializations and defaults
   ! tdstep of the predictor is 0.5 * tdstep_magnus
   ! F1a and F1b are used to extrapolate F3, then F3 is used to propagate rho.
   ! Afterwards, rho4 is copied into Pmat_vec in order to calculate F5.
   tdstep1 = timestep * 0.50D0
   F3      = (7.D0 / 4.D0) * F1b - (3.D0 / 4.D0) * F1a
   rho2t   = rho2

   call magnus(F3(:,:,1), rho2(:,:,1), rho4(:,:,1), M_in, NBCH, tdstep1, &
               factorial)
   if (OPEN) then
      call magnus(F3(:,:,2), rho2(:,:,2), rho4(:,:,2), M_in, NBCH, tdstep1, &
                  factorial)
      rho2t(:,:,1) = basechange(M_in, Xmat, rho4(:,:,1), Xtrans)
      rho2t(:,:,2) = basechange(M_in, Xmat, rho4(:,:,2), Xtrans)
      call sprepack_ctr('L', M, rhoalpha, rho2t(MTB+1:MTB+M,MTB+1:MTB+M,1))
      call sprepack_ctr('L', M, rhobeta , rho2t(MTB+1:MTB+M,MTB+1:MTB+M,2))
      Pmat_vec = rhoalpha + rhobeta
   else
      rho2t(:,:,1) = basechange(M_in, Xmat, rho4(:,:,1), Xtrans)
      call sprepack_ctr('L', M, Pmat_vec, rho2t(MTB+1:MTB+M,MTB+1:MTB+M,1))
   end if

   call int3lu(E2, Pmat_vec, Fmat_vec2, Fmat_vec, Gmat_vec, Ginv_vec, &
               Hmat_vec, open, MEMO)
   call g2g_solve_groups(0, Ex, 0)
   call field_calc(E1, time, Pmat_vec, Fmat_vec2, Fmat_vec, r, d, &
                   Iz, natom, ntatom, open)
   FBA = FON
   call spunpack('L', M, Fmat_vec, FBA(MTB+1:MTB+M,MTB+1:MTB+M,1))

   call fockbias_apply(time, FBA(MTB+1:MTB+M,MTB+1:MTB+M,1))
   FON(:,:,1) = basechange(M_in, Xtrans, FBA(:,:,1), Xmat)

   if (OPEN) then
      call spunpack('L', M, Fmat_vec2, FBA(MTB+1:MTB+M,MTB+1:MTB+M,2))
      call fockbias_apply(time, FBA(MTB+1:MTB+M,MTB+1:MTB+M,2))
      FON(:,:,2) = basechange(M_in, Xtrans, FBA(:,:,2), Xmat)
   end if

   deallocate(rho4, rho2t, F3, FBA)
end subroutine predictor

#ifdef CUBLAS
subroutine cumagnusfac(Fock, RhoOld, RhoNew, M, N, dt, factorial)
   ! Propagator based in Baker-Campbell-Hausdorff Nth-order formula.
   ! Everything is calculated in the ortonormal basis.
   ! Input : Fock(t+(deltat/2)), rho(t)
   ! Output: rho6 = rho(t+deltat)
   implicit none
   integer         , intent(in)  :: M, N
   double precision, intent(in)  :: Fock(M,M), factorial(N), dt
#ifdef TD_SIMPLE
   complex         , intent(in)  :: RhoOld(M,M)
   complex         , intent(out) :: RhoNew(M,M)
#else
   double complex  , intent(in)  :: RhoOld(M,M)
   double complex  , intent(out) :: RhoNew(M,M)
#endif
   external :: CUBLAS_INIT, CUBLAS_SHUTDOWN, CUBLAS_SET_MATRIX, &
               CUBLAS_GET_MATRIX, CUBLAS_ALLOC, CUBLAS_CAXPY,   &
               CUBLAS_CGEMM, CUBLAS_ZCOPY
   integer  :: CUBLAS_INIT, CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_ZCOPY,  &
               CUBLAS_GET_MATRIX, CUBLAS_ZGEMM, CUBLAS_ZAXPY, CUBLAS_CAXPY, &
               CUBLAS_CGEMM
   integer          :: stat, icount, jcount
   integer*8        :: devPOmega, devPPrev, devPNext, devPRho, devPScratch
   double precision :: Fact
#ifdef TD_SIMPLE
   complex          :: alpha , beta
   complex, allocatable :: Omega1(:,:)
   complex, parameter   :: ICMPLX = CMPLX(0.0D0, 1.0D0)
   integer, parameter   :: SIZEOF_COMPLEX = 8
#else
   double complex   :: alpha , beta
   double complex, allocatable :: Omega1(:,:)
   double complex, parameter   :: ICMPLX = DCMPLX(0.0D0, 1.0D0)
   integer       , parameter   :: SIZEOF_COMPLEX = 16
#endif
   allocate(Omega1(M,M))
   do icount = 1, M
   do jcount = 1, M
      Omega1(icount,jcount) = -1 * ICMPLX * Fock(icount,jcount) * dt
   enddo
   enddo
   stat = CUBLAS_ALLOC(M*M, SIZEOF_COMPLEX, devPOmega)
   stat = CUBLAS_ALLOC(M*M, SIZEOF_COMPLEX, devPPrev)
   stat = CUBLAS_ALLOC(M*M, SIZEOF_COMPLEX, devPNext)
   stat = CUBLAS_ALLOC(M*M, SIZEOF_COMPLEX, devPRho)

   stat = CUBLAS_SET_MATRIX(M, M, SIZEOF_COMPLEX, Omega1, M, devPOmega, M)
   if (stat .ne. 0) then
      write(*,'(A)') "  ERROR: Omega1 allocation failed (magnus_cublas)."
      call CUBLAS_SHUTDOWN()
      stop
   endif

   stat = CUBLAS_SET_MATRIX(M, M, SIZEOF_COMPLEX, RhoOld, M, devPPrev, M)
   if (stat .ne. 0) then
      write(*,'(A)') "  ERROR: RhoPrev allocation failed (magnus_cublas)."
      call CUBLAS_SHUTDOWN()
      stop
   endif

   stat = CUBLAS_SET_MATRIX(M, M, SIZEOF_COMPLEX, RhoOld, M, devPRho, M)
   if (stat .ne. 0) then
      write(*,'(A)') "  ERROR: RhoOld allocation failed (magnus_cublas)."
      call CUBLAS_SHUTDOWN()
      stop
   endif

   ! Density matrix propagation
   Rhonew = RhoOld
   alpha  = (1.0D0, 0.0D0)
   do icount = 1, N
      Fact  = factorial(icount)
      beta  = (0.0D0, 0.0D0)
#ifdef TD_SIMPLE
      alpha = CMPLX(Fact, 0.0D0)
      stat  = CUBLAS_CGEMM('N', 'N', M, M, M, alpha, devPOmega, M, devPPrev, &
                           M, beta, devPNext, M)
#else
      alpha = DCMPLX(Fact, 0.0D0)
      stat  = CUBLAS_ZGEMM('N', 'N', M, M, M, alpha, devPOmega, M, devPPrev, &
                           M, beta, devPNext, M)
#endif

      if (stat .ne. 0) then
         write(*,'(A)') " ERROR: ZGEM1 failed (magnus_cublas)"
         call magnus_shutdown(devPOmega, devPRho, devPNext, devPPrev)
         stop
      endif

      beta  = (-1.0D0, 0.0D0)
#ifdef TD_SIMPLE
      stat  = CUBLAS_CGEMM('N', 'N', M, M, M, alpha, devPPrev, M, devPOmega, &
                           M, beta, devPNext, M)
#else
      stat  = CUBLAS_ZGEMM('N', 'N', M, M, M, alpha, devPPrev, M, devPOmega, &
                           M, beta, devPNext, M)
#endif
      if (stat .ne. 0) then
         write(*,'(A)') " ERROR: ZGEM2 failed (magnus_cublas)"
         call magnus_shutdown(devPOmega, devPRho, devPNext, devPPrev)
         stop
      endif

#ifdef TD_SIMPLE
      stat = CUBLAS_CAXPY(M*M, CMPLX(1.0D0,0.0D0), devPNext, 1, devPRho, 1)
#else
      stat = CUBLAS_ZAXPY(M*M, DCMPLX(1.0D0,0.0D0), devPNext, 1, devPRho, 1)
#endif
      if (stat .ne. 0) then
         write(*,'(A)') " ERROR: CAXPY failed (magnus_cublas)"
         call magnus_shutdown(devPOmega, devPRho, devPNext, devPPrev)
         stop
      endif

      devPScratch = devPPrev
      devPPrev    = devPNext
      devPNext    = devPScratch
   enddo

   stat = CUBLAS_GET_MATRIX(M, M, SIZEOF_COMPLEX, devPRho, M, Rhonew, M)
   if (stat .ne. 0) then
      write(*,'(A)') " ERROR: Get_matrix failed (magnus_cublas)"
      call magnus_shutdown(devPOmega, devPRho, devPNext, devPPrev)
      stop
   endif

   call CUBLAS_FREE ( devPOmega )
   call CUBLAS_FREE ( devPRho )
   call CUBLAS_FREE ( devPNext )
   call CUBLAS_FREE ( devPPrev )
   deallocate(Omega1)
end subroutine cumagnusfac

subroutine cupredictor(F1a, F1b, FON, rho2, devPtrX, factorial, devPtrXc, &
                       timestep, time, M_in, MTB, dim3)
   ! Predictor-Corrector Cheng, V.Vooris.PhysRevB.2006.74.155112
   ! This routine receives: F1a, F1b, rho2
   ! And gives: F5 = F(t+(deltat/2))
   use garcha_mod   , only: NBCH, rhoalpha, rhobeta, OPEN, r, d, natom,      &
                            ntatom, Iz, MEMO, Fmat_vec, Fmat_vec2, Ginv_vec, &
                            Hmat_vec, Gmat_vec, Pmat_vec
   use field_data   , only: field
   use field_subs   , only: field_calc
   use cublasmath   , only: basechange_cublas
   use faint_cpu    , only: int3lu
   use fockbias_subs, only: fockbias_apply
   use basis_data   , only: M, Md

   implicit none
   ! MTB is only used in DFTB, it equals 0 otherwise.
   integer         , intent(in)    :: M_in, dim3, MTB
   integer*8       , intent(in)    :: devPtrX, devPtrXc
   double precision, intent(in)    :: timestep, factorial(NBCH), time
   double precision, intent(inout) :: F1a(M_in,M_in,dim3), F1b(M_in,M_in,dim3),&
                                      FON(M_in,M_in,dim3)
#ifdef TD_SIMPLE
   complex, intent(in)  :: rho2(M_in,M_in,dim3)
   complex, allocatable :: rho4(:,:,:), rho2t(:,:,:)
#else
   double complex, intent(in)  :: rho2(M_in,M_in,dim3)
   double complex, allocatable :: rho4(:,:,:), rho2t(:,:,:)
#endif
   integer :: M2
   double precision :: E2, tdstep1, Ex, E1
   double precision, allocatable :: F3(:,:,:), FBA(:,:,:)

   allocate(rho4(M_in,M_in,dim3), rho2t(M_in,M_in,dim3), F3(M_in,M_in,dim3), &
            FBA(M_in,M_in,dim3))

   M2 = 2 * M

   ! Initializations and defaults
   ! tdstep of the predictor is 0.5 * tdstep_magnus
   ! F1a and F1b are used to extrapolate F3, then F3 is used to propagate rho.
   ! Afterwards, rho4 is copied into Pmat_vec in order to calculate F5.
   tdstep1 = timestep * 0.50D0
   F3      = (7.D0 / 4.D0) * F1b - (3.D0 / 4.D0) * F1a
   rho2t   = rho2

   call cumagnusfac(F3(:,:,1), rho2(:,:,1), rho4(:,:,1), M_in, NBCH, tdstep1, &
                    factorial)
   if (OPEN) then
      call cumagnusfac(F3(:,:,2), rho2(:,:,2), rho4(:,:,2), M_in, NBCH,tdstep1,&
                       factorial)
      rho2t(:,:,1) = basechange_cublas(M_in, rho4(:,:,1), devPtrXc, 'inv')
      rho2t(:,:,2) = basechange_cublas(M_in, rho4(:,:,2), devPtrXc, 'inv')
      call sprepack_ctr('L', M, rhoalpha, rho2t(MTB+1:MTB+M,MTB+1:MTB+M,1))
      call sprepack_ctr('L', M, rhobeta , rho2t(MTB+1:MTB+M,MTB+1:MTB+M,2))
      Pmat_vec = rhoalpha + rhobeta
   else
      rho2t(:,:,1) = basechange_cublas(M_in, rho4(:,:,1), devPtrXc, 'inv')
      call sprepack_ctr('L', M, Pmat_vec, rho2t(MTB+1:MTB+M,MTB+1:MTB+M,1))
   end if

   call int3lu(E2, Pmat_vec, Fmat_vec2, Fmat_vec, Gmat_vec, Ginv_vec, &
               Hmat_vec, open, MEMO)
   call g2g_solve_groups(0, Ex, 0)
   call field_calc(E1, time, Pmat_vec, Fmat_vec2, Fmat_vec, r, d, &
                   Iz, natom, ntatom, open)
   FBA = FON
   call spunpack('L', M, Fmat_vec, FBA(MTB+1:MTB+M,MTB+1:MTB+M,1))

   call fockbias_apply(time, FBA(MTB+1:MTB+M,MTB+1:MTB+M,1))
   FON(:,:,1) = basechange_cublas(M_in, FBA(:,:,1), devPtrX, 'dir')

   if (OPEN) then
      call spunpack('L', M, Fmat_vec2, FBA(MTB+1:MTB+M,MTB+1:MTB+M,2))
      call fockbias_apply(time, FBA(MTB+1:MTB+M,MTB+1:MTB+M,2))
      FON(:,:,2) = basechange_cublas(M_in, FBA(:,:,2), devPtrX, 'dir')
   end if

   deallocate(rho4, rho2t, F3, FBA)
end subroutine cupredictor

subroutine magnus_cublas(Fock, RhoOld, RhoNew, M, N, dt, factorial)
   ! Propagator based in Baker-Campbell-Hausdorff Nth-order formula.
   ! Everything is calculated in the ortonormal basis.
   ! Input : Fock(t+(deltat/2)), rho(t)
   ! Output: rho6 = rho(t+deltat)
   implicit none
   integer         , intent(in)  :: M, N
   double precision, intent(in)  :: Fock(M,M), factorial(N), dt
   double complex  , intent(in)  :: RhoOld(M,M)
   double complex  , intent(out) :: RhoNew(M,M)

   external :: CUBLAS_INIT, CUBLAS_SHUTDOWN, CUBLAS_SET_MATRIX, &
               CUBLAS_GET_MATRIX, CUBLAS_ALLOC, CUBLAS_CAXPY,   &
               CUBLAS_CGEMM, CUBLAS_ZCOPY
   integer  :: CUBLAS_INIT, CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_ZCOPY,  &
               CUBLAS_GET_MATRIX, CUBLAS_ZGEMM, CUBLAS_ZAXPY, CUBLAS_CAXPY, &
               CUBLAS_CGEMM

   integer          :: stat, icount, jcount
   integer*8        :: devPOmega, devPPrev, devPNext, devPRho, devPScratch
   double precision :: Fact
#ifdef TD_SIMPLE
   complex   :: alpha , beta
   complex, allocatable :: Omega1(:,:)
   complex, parameter   :: ICMPLX = CMPLX(0.0D0, 1.0D0)
   integer, parameter   :: SIZEOF_COMPLEX = 8
#else
   double complex   :: alpha , beta
   double complex, allocatable :: Omega1(:,:)
   double complex, parameter   :: ICMPLX = DCMPLX(0.0D0, 1.0D0)
   integer       , parameter   :: SIZEOF_COMPLEX = 16
#endif
   allocate(Omega1(M,M))
   do icount = 1, M
   do jcount = 1, M
      Omega1(icount,jcount) = ICMPLX * Fock(icount,jcount) * dt
   enddo
   enddo
   stat = CUBLAS_ALLOC(M*M, SIZEOF_COMPLEX, devPOmega)
   stat = CUBLAS_ALLOC(M*M, SIZEOF_COMPLEX, devPPrev)
   stat = CUBLAS_ALLOC(M*M, SIZEOF_COMPLEX, devPNext)
   stat = CUBLAS_ALLOC(M*M, SIZEOF_COMPLEX, devPRho)

   stat = CUBLAS_SET_MATRIX(M, M, SIZEOF_COMPLEX, Omega1, M, devPOmega, M)
   if (stat .ne. 0) then
      write(*,'(A)') "  ERROR: Omega1 allocation failed (magnus_cublas)."
      call CUBLAS_SHUTDOWN()
      stop
   endif

   stat = CUBLAS_SET_MATRIX(M, M, SIZEOF_COMPLEX, RhoOld, M, devPPrev, M)
   if (stat .ne. 0) then
      write(*,'(A)') "  ERROR: RhoOld allocation failed (magnus_cublas)."
      call CUBLAS_SHUTDOWN()
      stop
   endif

   stat = CUBLAS_ZCOPY(M*M, devPPrev, 1,devPNext , 1)
   if (stat .ne. 0) then
      write(*,'(A)') "  ERROR: Matrix copy failed (magnus_cublas)."
      call CUBLAS_SHUTDOWN()
      stop
   endif

   ! Density matrix propagation
   Rhonew = RhoOld
   alpha  = (1.0D0, 0.0D0)
   do icount = 1, N
      Fact  = factorial(icount)
      beta   = (0.0D0, 0.0D0)
#ifdef TD_SIMPLE
      alpha = CMPLX(Fact, 0.0D0)
      stat  = CUBLAS_CGEMM('N', 'N', M, M, M, alpha, devPPrev, M, devPOmega,&
                           M, beta, devPNext, M)
#else
      alpha = DCMPLX(Fact, 0.0D0)
      stat  = CUBLAS_ZGEMM('N', 'N', M, M, M, alpha, devPPrev, M, devPOmega,&
                           M, beta, devPNext, M)
#endif
      if (stat .ne. 0) then
         write(*,'(A)') " ERROR: ZGEM1 failed (magnus_cublas)"
         call magnus_shutdown(devPOmega, devPRho, devPNext, devPPrev)
         stop
      endif

      beta   = (-1.0D0, 0.0D0)
#ifdef TD_SIMPLE
      alpha = CMPLX(Fact, 0.0D0)
      stat  = CUBLAS_CGEMM('N', 'N', M, M, M, alpha, devPOmega, M, devPPrev, &
                           M, beta, devPNext, M)
#else
      alpha = DCMPLX(Fact, 0.0D0)
      stat  = CUBLAS_ZGEMM('N', 'N', M, M, M, alpha, devPOmega, M, devPPrev, &
                           M, beta, devPNext, M)
#endif
      if (stat .ne. 0) then
         write(*,'(A)') " ERROR: ZGEM2 failed (magnus_cublas)"
         call magnus_shutdown(devPOmega, devPRho, devPNext, devPPrev)
         stop
      endif

      if (icount .ne. N) then
#ifdef TD_SIMPLE
         stat = CUBLAS_CAXPY(M*M, CMPLX(1.0D0,0.0D0), devPNext, 1, devPRho, 1)
#else
         stat = CUBLAS_ZAXPY(M*M, DCMPLX(1.0D0,0.0D0), devPNext, 1, devPRho, 1)
#endif
         if (stat .ne. 0) then
            write(*,'(A)') " ERROR: ZCOPY failed (magnus_cublas)"
            call magnus_shutdown(devPOmega, devPRho, devPNext, devPPrev)
            stop
         endif
      endif
      devPScratch = devPPrev
      devPPrev    = devPNext
      devPNext    = devPScratch
   enddo

   stat = CUBLAS_GET_MATRIX(M, M, SIZEOF_COMPLEX, devPRho, M, Rhonew, M)
   if (stat .ne. 0) then
      write(*,'(A)') " ERROR: Get_matrix failed (magnus_cublas)"
      call magnus_shutdown(devPOmega, devPRho, devPNext, devPPrev)
      stop
   endif

   call CUBLAS_FREE ( devPOmega )
   call CUBLAS_FREE ( devPRho )
   call CUBLAS_FREE ( devPNext )
   call CUBLAS_FREE ( devPPrev )
   deallocate(Omega1)
end subroutine magnus_cublas

subroutine magnus_shutdown(PTR1, PTR2, PTR3, PTR4)
   implicit none
   external :: CUBLAS_FREE, CUBLAS_SHUTDOWN
   integer*8, intent(in) :: PTR1, PTR2, PTR3, PTR4

   call CUBLAS_FREE(PTR1)
   call CUBLAS_FREE(PTR2)
   call CUBLAS_FREE(PTR3)
   call CUBLAS_FREE(PTR4)
   call CUBLAS_SHUTDOWN()

end subroutine magnus_shutdown
#endif
end module propagators
