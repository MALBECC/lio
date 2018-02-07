!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% MAGNUS.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains the routine for an N-order magnus propagation.            !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine magnus(fock, rhoOld, rhoNew, M, N, dt, factorial)
    ! Input:  Fock(t+(deltat/2)), rho(t)
    ! Output: rho6 = rho(t+deltat)
    implicit none
    integer  , intent(in)  :: M, N
    real*8   , intent(in)  :: Fock(M,M), dt, factorial(N)
#ifdef TD_SIMPLE
    complex*8 , intent(in)  :: RhoOld(M,M)
    complex*8 , intent(out) :: RhoNew(M,M)
    complex*8 , parameter   :: icmplx=CMPLX(0.0D0,1.0D0)
    complex*8 , allocatable :: Scratch(:,:), ConmPrev(:,:), ConmNext(:,:), &
                               Omega1(:,:)
#else
    complex*16, intent(in)  :: RhoOld(M,M)
    complex*16, intent(out) :: RhoNew(M,M)
    complex*16, parameter   :: icmplx=CMPLX(0.0D0,1.0D0)
    complex*16, allocatable :: Scratch(:,:), ConmPrev(:,:), ConmNext(:,:), &
                               Omega1(:,:)
#endif
    integer                :: ii

    ! Variable initializations.
    allocate(ConmPrev(M,M), ConmNext(M,M), Omega1(M,M), Scratch(M,M))
    ConmPrev = RhoOld
    RhoNew   = RhoOld

    ! Omega1 (=W) instantiation.
    Omega1 = (-1)*(icmplx)*(Fock)*(dt)

    ! Density matrix propagation
    do ii=1, N
        ConmNext = MATMUL(Omega1, ConmPrev)
        Scratch  = MATMUL(ConmPrev, Omega1)
        ConmNext = ConmNext - Scratch
        RhoNew   = RhoNew + factorial(ii)*ConmNext
        ConmPrev = ConmNext
    enddo

    deallocate(ConmPrev, ConmNext, Omega1, Scratch)
    return
end subroutine magnus
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
