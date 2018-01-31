!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine predictor(F1a, F1b, FON, rho2, factorial, Xmat,
     >                      Xtrans, fxx, fyy, fzz,g)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This routine recives: F1a,F1b,rho2
! And gives: F5 = F(t+(deltat/2))
       use garcha_mod , only: M, RMM, NBCH, field
       use td_data    , only: tdstep
       use mathsubs   , only: basechange
       use faint_cpu77, only: int3lu, intfld
       implicit none
       REAL*8,intent(inout) :: F1a(M,M),F1b(M,M),FON(M,M), Xmat(M,M)
       REAL*8,intent(in) :: Xtrans(M,M), fxx, fyy, fzz,g
       REAL*8, intent(in) :: factorial(NBCH)
       REAL*8,allocatable :: F3(:,:),FBA(:,:)
#ifdef TD_SIMPLE
       COMPLEX*8, intent(in) :: rho2(M,M)
       complex*8,allocatable :: rho4(:,:),rho2t(:,:)
#else
      COMPLEX*16, intent(in) :: rho2(M,M)
      complex*16,allocatable :: rho4(:,:),rho2t(:,:)
#endif
       integer :: i,j,k,kk, M2, M5, MM
       real*8 :: E2, tdstep1, Ex
!------------------------------------------------------------------------------!
       ALLOCATE(rho4(M,M),rho2t(M,M),F3(M,M),FBA(M,M))
c
       M2=2*M
       MM=M*(M+1)/2
c now S, F also uses the same position after S was used
       M5=1 + 2*MM

c Initializations/Defaults
! tdstep predictor is 0.5*tdstep magnus
       tdstep1=tdstep*0.50D0
! Step 1: Matrix F1a and F1b are used to extrapolate F3
       F3=(7.D0/4.D0)*F1b-(3.D0/4.D0)*F1a
! Step2: F3 is used to propagate rho2 to rho4
       rho2t=rho2
       call magnus(F3,rho2,rho4,M,NBCH,tdstep1,factorial)
       rho2t = basechange(M,Xmat,rho4,Xtrans)
! Paso3: Escribimos rho4 en el RMM para poder obtener F5 en el siguiente paso.
! Step3: rho4 is copied to RMM(1,2,3,...,MM)
      call sprepack_ctr('L',M,RMM,rho2t)
! Step4: Density matrix 4 is used to calculate F5
       call int3lu(E2)
       call g2g_solve_groups(0,Ex,0)
       if (field) call intfld(g,Fxx,Fyy,Fzz)
       call spunpack('L',M,RMM(M5),FBA)
       FON=basechange(M,Xtrans,FBA,Xmat)
       DEALLOCATE(rho4,rho2t,F3,FBA)
       RETURN;END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
