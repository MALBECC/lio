!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine cupredictor_op_DZ(F1a_a,F1b_a,F1a_b,F1b_b,FON_a,FON_b,
     > rho2_a,rho2_b,factorial,devPtrX,devPtrXc, time)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This routine recives: F1a,F1b,rho2
! And gives: F5 = F(t+(deltat/2))
       use garcha_mod
       use field_data, only: field
       use field_subs, only: field_calc
       use faint_cpu77, only: int3lu
       REAL*8,intent(inout) :: F1a_a(M,M),F1b_a(M,M),
     > F1a_b(M,M),F1b_b(M,M),FON_a(M,M),FON_b(M,M)
       integer*8,intent(in) :: devPtrX,devPtrXc
       REAL*8, intent(in) :: factorial(NBCH), time
       REAL*8,allocatable :: F3(:,:),FBA(:,:)
       integer :: i,j,k,kk,stat
       real*8 :: E2, tdstep1, E1
       COMPLEX*16, intent(in) :: rho2_a(M,M),rho2_b(M,M)
       COMPLEX*16,allocatable :: rho4(:,:),rho2t(:,:)
       integer*8 devPtrScratch1
!------------------------------------------------------------------------------!
       ALLOCATE(rho4(M,M),rho2t(M,M),F3(M,M),FBA(M,M))
c
       M2=2*M
       MM=M*(M+1)/2
c first i
       M1=1
c now Fold
       M3=M1+MM
c now S, F also uses the same position after S was used
       M5=M3+MM
c now G
! tdstep predictor is 0.5*tdstep magnus
       tdstep1=tdstep*0.50D0
!
       F3=(7.D0/4.D0)*F1b_a-(3.D0/4.D0)*F1a_a
!
!       call magnus(F3,rho2_a,rho2t,M,NBCH,tdstep1,factorial)
       call cumagnusfac(F3,rho2_a,rho2t,M,NBCH,tdstep1,factorial)
!
!       call matmulnanoc(rho2t,xtrans,rho4,M)
!       call rho_transform(rho2t,devPtrX,rho4,M)
!       call complex_rho_on_to_ao(rho2t,devPtrXc,rho4,M)
       rho4=basechange_cublas(M,rho2t,devPtrXc,'inv')
       call sprepack_ctr('L',M,rhoalpha,rho4)
       call sprepack_ctr('L',M,RMM,rho4)
!
       F3=(7.D0/4.D0)*F1b_b-(3.D0/4.D0)*F1a_b
!
!       call magnus(F3,rho2_b,rho2t,M,NBCH,tdstep1,factorial)
        call cumagnusfac(F3,rho2_b,rho2t,M,NBCH,tdstep1,factorial)
!
       rho4=0
!       call matmulnanoc(rho2t,xtrans,rho4,M)
!       call rho_transform(rho2t,devPtrX,rho4,M)
!       call complex_rho_on_to_ao(rho2t,devPtrXc,rho4,M)
       rho4=basechange_cublas(M,rho2t,devPtrXc,'inv')
       call sprepack_ctr('L',M,rhobeta,rho4)
       DO i=1,MM
          RMM(i)=RMM(i)+rhobeta(i)
       ENDDO
! Step4: Density matrix 4 is used to calculate F5
       call int3lu(E2)
       call g2g_solve_groups(0,Ex,0)
       call field_calc(E1, time)
!
       call spunpack('L',M,RMM(M5),FBA)
!       call matmulnano(FBA,X,FON_a,M)
!       call fock_ao_to_on(FBA,devPtrX,FON_a,M)
       FON_a=basechange_cublas(M,FBA,devPtrX,'dir')
       call spunpack('L',M,RMM(M3),FBA)
!       call matmulnano(FBA,X,FON_b,M)
!       call fock_ao_to_on(FBA,devPtrX,FON_b,M)
       FON_b=basechange_cublas(M,FBA,devPtrX,'dir')
       DEALLOCATE(rho4,rho2t,F3,FBA)
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine cupredictor_op_DC(F1a_a,F1b_a,F1a_b,F1b_b,FON_a,FON_b,
     > rho2_a,rho2_b,factorial,devPtrX,devPtrXc, time)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This routine recives: F1a,F1b,rho2
! And gives: F5 = F(t+(deltat/2))
       use garcha_mod
       use field_data, only: field
       use field_subs, only: field_calc
       REAL*8,intent(inout) :: F1a_a(M,M),F1b_a(M,M),
     > F1a_b(M,M),F1b_b(M,M),FON_a(M,M),FON_b(M,M)
       integer*8,intent(in) :: devPtrX,devPtrXc
       REAL*8, intent(in) :: factorial(NBCH), time
       REAL*8,allocatable :: F3(:,:),FBA(:,:)
       integer :: i,j,k,kk,stat
       real*8 :: E2, tdstep1, E1
       COMPLEX*8, intent(in) :: rho2_a(M,M),rho2_b(M,M)
       COMPLEX*8,allocatable :: rho4(:,:),rho2t(:,:)
       integer*8 devPtrScratch1
!------------------------------------------------------------------------------!
       ALLOCATE(rho4(M,M),rho2t(M,M),F3(M,M),FBA(M,M))
c
       M2=2*M
       MM=M*(M+1)/2
c first i
       M1=1
c now Fold
       M3=M1+MM
c now S, F also uses the same position after S was used
       M5=M3+MM
c now G
! tdstep predictor is 0.5*tdstep magnus
       tdstep1=tdstep*0.50D0
!
       F3=(7.D0/4.D0)*F1b_a-(3.D0/4.D0)*F1a_a
!
!       call magnus(F3,rho2_a,rho2t,M,NBCH,tdstep1,factorial)
       call cumagnusfac(F3,rho2_a,rho2t,M,NBCH,tdstep1,factorial)
!
!       call matmulnanoc(rho2t,xtrans,rho4,M)
!       call rho_transform(rho2t,devPtrX,rho4,M)
!       call complex_rho_on_to_ao(rho2t,devPtrXc,rho4,M)
       rho4=basechange_cublas(M,rho2t,devPtrXc,'inv')
       call sprepack_ctr('L',M,rhoalpha,rho4)
       call sprepack_ctr('L',M,RMM,rho4)
!
       F3=(7.D0/4.D0)*F1b_b-(3.D0/4.D0)*F1a_b
!
!       call magnus(F3,rho2_b,rho2t,M,NBCH,tdstep1,factorial)
        call cumagnusfac(F3,rho2_b,rho2t,M,NBCH,tdstep1,factorial)
!
       rho4=0
!       call matmulnanoc(rho2t,xtrans,rho4,M)
!       call rho_transform(rho2t,devPtrX,rho4,M)
!       call complex_rho_on_to_ao(rho2t,devPtrXc,rho4,M)
       rho4=basechange_cublas(M,rho2t,devPtrXc,'inv')
       call sprepack_ctr('L',M,rhobeta,rho4)
       DO i=1,MM
          RMM(i)=RMM(i)+rhobeta(i)
       ENDDO
! Step4: Density matrix 4 is used to calculate F5
       call int3lu(E2)
       call g2g_solve_groups(0,Ex,0)
       call field_calc(E1, time)
!
       call spunpack('L',M,RMM(M5),FBA)
!       call matmulnano(FBA,X,FON_a,M)
!       call fock_ao_to_on(FBA,devPtrX,FON_a,M)
       FON_a=basechange_cublas(M,FBA,devPtrX,'dir')
       call spunpack('L',M,RMM(M3),FBA)
!       call matmulnano(FBA,X,FON_b,M)
!       call fock_ao_to_on(FBA,devPtrX,FON_b,M)
       FON_b=basechange_cublas(M,FBA,devPtrX,'dir')
       DEALLOCATE(rho4,rho2t,F3,FBA)
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
