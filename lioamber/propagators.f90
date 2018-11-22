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
    complex, parameter   :: icmplx = CMPLX(0.0D0,1.0D0)
    complex, allocatable :: Scratch(:,:), ConmPrev(:,:), ConmNext(:,:), &
                            Omega1(:,:)
#else
    double complex  , intent(in)  :: RhoOld(M,M)
    double complex  , intent(out) :: RhoNew(M,M)
    double complex, parameter   :: icmplx = CMPLX(0.0D0,1.0D0)
    double complex, allocatable :: Scratch(:,:), ConmPrev(:,:), ConmNext(:,:), &
                                   Omega1(:,:)
#endif
    integer :: icount

    ! Variable initializations.
    allocate(ConmPrev(M,M), ConmNext(M,M), Omega1(M,M), Scratch(M,M))
    ConmPrev = RhoOld
    RhoNew   = RhoOld

    ! Omega1 (=W) instantiation.
    Omega1 = -1 * icmplx * Fock * dt

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
   ! This routine recives: F1a,F1b,rho2
   ! And gives: F5 = F(t+(deltat/2))
   use garcha_mod   , only: M, RMM, NBCH, rhoalpha, rhobeta, OPEN, Md, r, d, &
                            natom, ntatom, Iz
   use field_data   , only: field
   use field_subs   , only: field_calc
   use mathsubs     , only: basechange
   use faint_cpu    , only: int3lu
   use fockbias_subs, only: fockbias_apply

   implicit none
   integer         , intent(in)    :: M_in, dim3
   integer         , intent(in)    :: MTB ! Only used in DFTB, =0 otherwise.
   double precision, intent(in)    :: Xtrans(M_in,M_in), timestep
   double precision, intent(in)    :: factorial(NBCH), time
   double precision, intent(inout) :: F1a(M_in,M_in,dim3), F1b(M_in,M_in,dim3),&
                                      FON(M_in,M_in,dim3), Xmat(M_in,M_in)
#ifdef TD_SIMPLE
   complex, intent(in)  :: rho2(M_in,M_in,dim3)
   complex, allocatable :: rho4(:,:,:), rho2t(:,:,:)
#else
   double complex, intent(in)  :: rho2(M_in,M_in,dim3)
   double complex, allocatable :: rho4(:,:,:), rho2t(:,:,:)
#endif
   integer :: i,j,k,kk, M2, M3, M5, MM, MMd, M11, M7, M9
   double precision :: E2, tdstep1, Ex, E1
   double precision, allocatable :: F3(:,:,:), FBA(:,:,:)

   allocate(rho4(M_in,M_in,dim3), rho2t(M_in,M_in,dim3), F3(M_in,M_in,dim3), &
            FBA(M_in,M_in,dim3))

   M2 = 2 * M    ;  MM  = M * (M +1)/2; MMd = Md * (Md+1)/2; M3  = MM +1
   M5 = 2 * MM +1;  M7  = M5 + MM     ; M9  = M7 + MMd     ; M11 = M9 + MMd

   ! Initializations and defaults
   ! tdstep of the predictor is 0.5 * tdstep_magnus
   ! F1a and F1b are used to extrapolate F3, then F3 is used to propagate rho.
   ! Afterwards, rho4 is copied into RMM in order to calculate F5.
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
      RMM(1:MM) = rhoalpha + rhobeta
   else
      rho2t(:,:,1) = basechange(M_in, Xmat, rho4(:,:,1), Xtrans)
      call sprepack_ctr('L', M, RMM, rho2t(MTB+1:MTB+M,MTB+1:MTB+M,1))
   end if

   call int3lu(E2, RMM(1:MM), RMM(M3:M3+MM), RMM(M5:M5+MM), RMM(M7:M7+MMd), &
               RMM(M9:M9+MMd), RMM(M11:M11+MMd), open)
   call g2g_solve_groups(0, Ex, 0)
   call field_calc(E1, time, RMM(M3:M3+MM), RMM(M5:M5+MM), r, d, Iz, natom, &
                   ntatom, open)
   FBA = FON
   call spunpack('L', M, RMM(M5), FBA(MTB+1:MTB+M,MTB+1:MTB+M,1))

   call fockbias_apply(time, FBA(MTB+1:MTB+M,MTB+1:MTB+M,1))
   FON(:,:,1) = basechange(M_in, Xtrans, FBA(:,:,1), Xmat)

   if (OPEN) then
      call spunpack('L', M, RMM(M3), FBA(MTB+1:MTB+M,MTB+1:MTB+M,2))
      call fockbias_apply(time, FBA(MTB+1:MTB+M,MTB+1:MTB+M,2))
      FON(:,:,2) = basechange(M_in, Xtrans, FBA(:,:,2), Xmat)
   end if

   deallocate(rho4, rho2t, F3, FBA)
end subroutine predictor
end module propagators
