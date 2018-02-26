!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine predictor(F1a, F1b, FON, rho2, factorial, Xmat,
     >                      Xtrans, timestep, time, M_in, MTB)
      ! This routine recives: F1a,F1b,rho2
      ! And gives: F5 = F(t+(deltat/2))
       use garcha_mod , only: M, RMM, NBCH
       use field_data , only: field
       use field_subs , only: field_calc
       use mathsubs   , only: basechange
       use faint_cpu77, only: int3lu
       implicit none
       REAL*8,intent(inout) :: F1a(M_in,M_in),F1b(M_in,M_in),
     >                         FON(M_in,M_in), Xmat(M_in,M_in)
       REAL*8,intent(in)  :: Xtrans(M_in,M_in), timestep
       REAL*8,intent(in)  :: factorial(NBCH), time
       REAL*8,allocatable :: F3(:,:), FBA(:,:)
#ifdef TD_SIMPLE
       COMPLEX*8, intent(in) :: rho2(M_in,M_in)
       complex*8,allocatable :: rho4(:,:),rho2t(:,:)
#else
      COMPLEX*16, intent(in) :: rho2(M_in,M_in)
      complex*16,allocatable :: rho4(:,:),rho2t(:,:)
#endif
       integer, intent(in)   :: M_in
       integer :: i,j,k,kk, M2, M5, MM
       real*8 :: E2, tdstep1, Ex, E1
       !DFTB: MTB variable is used for DFTB calculations otherwise equal to 0
       integer, intent(in)   :: MTB
!------------------------------------------------------------------------------!
       ALLOCATE(rho4(M_in,M_in),rho2t(M_in,M_in),F3(M_in,M_in),
     >          FBA(M_in,M_in))
c
       M2 = 2*M
       MM = M*(M+1)/2
       M5 = 1 + 2*MM

c Initializations/Defaults
! tdstep predictor is 0.5*tdstep magnus
       tdstep1=timestep*0.50D0
! Step 1: Matrix F1a and F1b are used to extrapolate F3
       F3=(7.D0/4.D0)*F1b-(3.D0/4.D0)*F1a
! Step2: F3 is used to propagate rho2 to rho4
       rho2t=rho2
       call magnus(F3,rho2,rho4,M_in,NBCH,tdstep1,factorial)
       rho2t = basechange(M_in,Xmat,rho4,Xtrans)
! Paso3: Escribimos rho4 en el RMM para poder obtener F5 en el siguiente paso.
! Step3: rho4 is copied to RMM(1,2,3,...,MM)
      call sprepack_ctr('L',M,RMM,rho2t(MTB+1:MTB+M,MTB+1:MTB+M))
! Step4: Density matrix 4 is used to calculate F5

       call int3lu(E2)
       call g2g_solve_groups(0,Ex,0)
       call field_calc(E1, time)
       FBA=FON
       call spunpack('L',M,RMM(M5),FBA(MTB+1:MTB+M,MTB+1:MTB+M))
       FON=basechange(M_in,Xtrans,FBA,Xmat)
       DEALLOCATE(rho4,rho2t,F3,FBA)
       RETURN;END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
