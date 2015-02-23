!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine predictor(F1a,F1b,FON,rho2,xtrans,factorial)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This routine recives: F1a,F1b,rho2
! And gives: F5 = F(t+(deltat/2))      
       use garcha_mod
       use mathsubs
       REAL*8,intent(inout) :: F1a(M,M),F1b(M,M),FON(M,M)
       REAL*8,intent(in) :: Xtrans(M,M)
       REAL*8, intent(in) :: factorial(NBCH)
       REAL*8,allocatable :: F3(:,:),FBA(:,:)
       COMPLEX*8, intent(in) :: rho2(M,M)
       complex*8,allocatable :: rho4(:,:),rho2t(:,:)
       integer :: i,j,k,kk
       real*8 :: E2, tdstep1
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
       M7=M5+MM
c now Gm
       M9=M7+MMd
c now H
       M11=M9+MMd
c W ( eigenvalues ), also this space is used in least squares
       M13=M11+MM
c aux ( vector for ESSl)
       M15=M13+M
c Least squares
       M17=M15+MM
c vectors of MO
       M18=M17+MMd
c weights (in case of using option )
       M19=M18+M*NCO
c
* RAM storage of two-electron integrals (if MEMO=T)
       M20 = M19 + natom*50*Nang
c Initializations/Defaults
! tdstep predictor is 0.5*tdstep magnus
       tdstep1=tdstep*0.50D0
! Step 1: Matrix F1a and F1b are used to extrapolate F3
       F3=(7.D0/4.D0)*F1b-(3.D0/4.D0)*F1a
! Step2: F3 is used to propagate rho2 to rho4
       rho2t=rho2
       call magnus(F3,rho2,rho4,M,NBCH,tdstep1,factorial)
! Paso3: Escribimos rho4 en el RMM para poder obtener F5 en el siguiente paso.
! Step3: rho4 is copied to RMM(1,2,3,...,MM)
       rho2t=basechange(M,X,rho4,Xtrans)
       do j=1,M
          do k=j,M
             if(j.eq.k) then
                RMM(k+(M2-j)*(j-1)/2)=REAL(rho2t(j,k))
             else
                RMM(k+(M2-j)*(j-1)/2)=REAL(rho2t(j,k))*2
             endif
          enddo
       enddo
! Step4: Density matrix 4 is used to calculate F5
       call int3lu(E2)
       call g2g_solve_groups(0,Ex,0)
       do j=1,M
          do k=1,j
             FBA(j,k)=RMM(M5+j+(M2-k)*(k-1)/2-1)
          enddo
          do k=j+1,M
             FBA(j,k)=RMM(M5+k+(M2-j)*(j-1)/2-1)
          enddo
       enddo
       FON=basechange(M,Xtrans,FBA,X)
       DEALLOCATE(rho4,rho2t,F3,FBA)
       RETURN;END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
