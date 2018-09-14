!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine predictor(F1a, F1b, FON, rho2, factorial, Xmat,
     >                      Xtrans, timestep, time, M_in, MTB, dim3)
      ! This routine recives: F1a,F1b,rho2
      ! And gives: F5 = F(t+(deltat/2))
       use garcha_mod , only: M, RMM, NBCH, rhoalpha, rhobeta, OPEN,
     >                        Md, r, d, Iz, natom, ntatom
       use field_data , only: field
       use field_subs , only: field_calc
       use mathsubs   , only: basechange
       use faint_cpu  , only: int3lu
       use fockbias_subs , only: fockbias_apply
       implicit none
       integer, intent(in)   :: M_in, dim3
       REAL*8,intent(inout) :: F1a(M_in,M_in,dim3),F1b(M_in,M_in,dim3),
     >                         FON(M_in,M_in,dim3), Xmat(M_in,M_in)
       REAL*8,intent(in)  :: Xtrans(M_in,M_in), timestep
       REAL*8,intent(in)  :: factorial(NBCH), time
       REAL*8,allocatable :: F3(:,:,:), FBA(:,:,:)
#ifdef TD_SIMPLE
       COMPLEX*8, intent(in) :: rho2(M_in,M_in,dim3)
       complex*8,allocatable :: rho4(:,:,:),rho2t(:,:,:)
#else
      COMPLEX*16, intent(in) :: rho2(M_in,M_in,dim3)
      complex*16,allocatable :: rho4(:,:,:),rho2t(:,:,:)
#endif
       integer :: i,j,k,kk, M2, M3, M5, MM, MMd, M11, M7, M9
       real*8 :: E2, tdstep1, Ex, E1
       !DFTB: MTB variable is used for DFTB calculations otherwise equal to 0
       integer, intent(in)   :: MTB
!------------------------------------------------------------------------------!
       ALLOCATE(rho4(M_in,M_in,dim3),rho2t(M_in,M_in,dim3),
     >          F3(M_in,M_in,dim3), FBA(M_in,M_in,dim3))
c
       M2 = 2*M
       MM = M*(M+1)/2
       MMd = Md*(Md+1)/2
       M3 = 1+MM
       M5 = 1 + 2*MM
       M7 = M5+MM ! G matrix
       M9 = M7+MMd ! G inverted
       M11= M9+MMd ! Hmat

c Initializations/Defaults
! tdstep predictor is 0.5*tdstep magnus
       tdstep1=timestep*0.50D0
! Step 1: Matrix F1a and F1b are used to extrapolate F3
       F3=(7.D0/4.D0)*F1b-(3.D0/4.D0)*F1a
! Step2: F3 is used to propagate rho2 to rho4
       rho2t=rho2
       call magnus(F3(:,:,1),rho2(:,:,1),rho4(:,:,1),M_in,NBCH,tdstep1,
     >             factorial)

       if (OPEN) then
         call magnus(F3(:,:,2),rho2(:,:,2),rho4(:,:,2),M_in,NBCH,
     >               tdstep1,factorial)
         rho2t(:,:,1) = basechange(M_in,Xmat,rho4(:,:,1),Xtrans)
         rho2t(:,:,2) = basechange(M_in,Xmat,rho4(:,:,2),Xtrans)
! Paso3open: Escribimos rho4 en rhoalpha y rhobeta para poder obtener
!            F5 en el siguiente paso.
           call sprepack_ctr('L',M,rhoalpha,
     >                       rho2t(MTB+1:MTB+M,MTB+1:MTB+M,1))
           call sprepack_ctr('L',M,rhobeta,
     >                       rho2t(MTB+1:MTB+M,MTB+1:MTB+M,2))
           RMM(1:MM) = rhoalpha + rhobeta
       else
           rho2t(:,:,1) = basechange(M_in,Xmat,rho4(:,:,1),Xtrans)
! Paso3: Escribimos rho4 en el RMM para poder obtener F5 en el siguiente paso.
! Step3: rho4 is copied to RMM(1,2,3,...,MM)
           call sprepack_ctr('L',M,RMM,
     >                       rho2t(MTB+1:MTB+M,MTB+1:MTB+M,1))
       end if
! Step4: Density matrix 4 is used to calculate F5
       call int3lu(E2, RMM(1:MM), RMM(M3:M3+MM), RMM(M5:M5+MM),
     >             RMM(M7:M7+MMd), RMM(M9:M9+MMd), RMM(M11:M11+MMd),
     >             open)
       call g2g_solve_groups(0,Ex,0)
       call field_calc(E1, time, RMM(M3:M3+MM), RMM(M5:M5+MM), r, d,
     > Iz, natom, ntatom, open)
       FBA=FON
       call spunpack('L',M,RMM(M5),FBA(MTB+1:MTB+M,MTB+1:MTB+M,1))

!Fockbias:
       call fockbias_apply(time, FBA(MTB+1:MTB+M,MTB+1:MTB+M,1))

       FON(:,:,1)=basechange(M_in,Xtrans,FBA(:,:,1),Xmat)

       if (OPEN) then
          call spunpack('L',M,RMM(M3),FBA(MTB+1:MTB+M,MTB+1:MTB+M,2))
!Fockbias:
          call fockbias_apply(time, FBA(MTB+1:MTB+M,MTB+1:MTB+M,2))

          FON(:,:,2)=basechange(M_in,Xtrans,FBA(:,:,2),Xmat)
       end if

       DEALLOCATE(rho4,rho2t,F3,FBA)
       RETURN;END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
