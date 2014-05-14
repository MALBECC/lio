!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE magnus(Fock,RhoOld,RhoNew,M,N,dt,factorial)
!------------------------------------------------------------------------------!
! Magnus propagator (N order)
! Entrada: Fock(t+(deltat/2)), rho(t)
! Salida:  rho6=rho(t+deltat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       IMPLICIT NONE
       INTEGER,INTENT(IN)                   :: M,N
       REAL*8,INTENT(IN)                    :: Fock(M,M),dt
       COMPLEX*8,INTENT(IN)                 :: RhoOld(M,M)
       COMPLEX*8,INTENT(OUT)                :: RhoNew(M,M)
       COMPLEX*8,ALLOCATABLE,DIMENSION(:,:) :: ConmPrev,ConmNext,Omega1
       COMPLEX*8,ALLOCATABLE,DIMENSION(:,:) :: Scratch
       INTEGER                              :: ii
       COMPLEX*8,PARAMETER                  :: icmplx=CMPLX(0.0D0,1.0D0)
       REAL*8, INTENT(IN)                   :: factorial(N)
!------------------------------------------------------------------------------!
! Variable initializations
       ALLOCATE(ConmPrev(M,M),ConmNext(M,M),Omega1(M,M))
       ALLOCATE(Scratch(M,M))
       ConmPrev=RhoOld
       RhoNew=RhoOld
! Omega1 (=W) instantiation
       Omega1=(-1)*(icmplx)*(Fock)*(dt)
! Density matrix propagation
       DO ii=1,N
         ConmNext=MATMUL(Omega1,ConmPrev)
         Scratch=MATMUL(ConmPrev,Omega1)
         ConmNext=ConmNext-Scratch
         RhoNew=RhoNew+factorial(ii)*ConmNext
         ConmPrev=ConmNext
       ENDDO
       DEALLOCATE(ConmPrev,ConmNext,Omega1,Scratch)
       RETURN;END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
