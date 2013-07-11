         subroutine predictor(F1a,F1b,FON,rho2,xtrans)
         use garcha_mod    

!Predictor-Corrector Cheng, V.Vooris.PhysRevB.2006.74.155112
! Esta rutina recibe: F1a,F1b,rho2
! Tira: F5 = F(t+(deltat/2))      

       REAL*8 , intent(inout) :: F1a(M,M),F1b(M,M),FON(M,M),Xtrans(M,M)
       REAL*8 :: F3(M,M),FBA(M,M)
       COMPLEX*8, intent(inout) :: rho2(M,M)
       complex*8 :: rho4(M,M)
       integer :: i,j,k,kk
       real*8 :: E2
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
c
c Initializations/Defaults
c xmm es la primer matriz de (M,M) en el 
! Paso1: Con las matrices pasadas F1a y F1b extrapolamos a F3----> Extrapolacion
       F3=(7.D0/4.D0)*F1b-((3.D0/4.D0)*F1a)
! Paso2: Usando H3, la matriz densidad rho2 es propagada a rho4----> Prediccion
       call magnus(F3,rho2,rho4,M,NBCH)
! Paso3: Escribimos rho4 en el RMM para poder obtener F5 en el siguiente paso.
        call matmulnanoc(rho4,xtrans,rho2,M)
        do j=1,M
         do k=j,M
          if(j.eq.k) then
           RMM(k+(M2-j)*(j-1)/2)=REAL(rho2(j,k))
          else
           RMM(k+(M2-j)*(j-1)/2)=REAL(rho2(j,k))*2
          endif
         enddo
        enddo
c        do i=1,MM
c         write(77,*) RMM(i)
c         enddo
c         stop
! Paso4: La matriz densidad 4 es usada para calcular F5------> Corrector
         call int3lu(E2)
         call g2g_solve_groups(0,Ex,0)
         write(*,*) 'E2=',E2,Ex         
!-------------------------Escritura de fock cuadrada--------------------------------------
!-----------Parte de abajo a la izquierda(incluyendo terminos diagonales)-----------------
         do j=1,M
          do k=1,j
!
           FBA(j,k)=RMM(M5+j+(M2-k)*(k-1)/2-1)
!
          enddo
!-----------Parte de arriba a la derecha de la matriz (sin incluir terminos diagonales)---
          do k=j+1,M
!
           FBA(j,k)=RMM(M5+k+(M2-j)*(j-1)/2-1)
!
          enddo
         enddo
!--------------Ahora tenemos F5 transformada en base de ON y en su forma cuadrada---------------
           call matmulnano(FBA,X,FON,M)

       return
       end


























