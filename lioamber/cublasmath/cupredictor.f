!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine cupredictor_DZ(F1a,F1b,FON,rho2,devPtrX,factorial,
     > devPtrXc, timestep, time, M_in, MTB)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!Predictor-Corrector Cheng, V.Vooris.PhysRevB.2006.74.155112
! Esta rutina recibe: F1a,F1b,rho2
! Tira: F5 = F(t+(deltat/2))
      use garcha_mod
      use field_data, only: field
      use field_subs, only: field_calc
      use mathsubs
      use general_module
      use faint_cpu77, only: int3lu
       implicit real*8 (a-h,o-z)
!       IMPLICIT REAL*8 (a-h,o-z)
       REAL*8,intent(inout) :: F1a(M_in,M_in),F1b(M_in,M_in)
       REAL*8,intent(inout) :: FON(M_in,M_in)
       integer*8,intent(in) :: devPtrX,devPtrXc
       integer, intent(in)    :: M_in
       integer, intent(in)    :: MTB
!       integer :: i,j,k,kk,stat,M1,M2,MM,M5,M7,M9,MMD,M11,M13,M15,M17
!       integer :: M19,M20,M3,M18
       integer :: i,j,k,kk,stat
       REAL*8,allocatable :: F3(:,:),FBA(:,:)
       real*8 :: E1, E2, tdstep1
      external CUBLAS_INIT, CUBLAS_SHUTDOWN
      integer CUBLAS_INIT
      REAL*8,intent(in) :: factorial(NBCH), timestep, time
       COMPLEX*16, intent(in) :: rho2(M_in,M_in)
       COMPLEX*16,allocatable :: rho4(:,:),rho2t(:,:)
!-----------------------------------------------------------------------------n
      ALLOCATE(rho4(M_in,M_in),rho2t(M_in,M_in),F3(M_in,M_in),
     >         FBA(M_in,M_in))
      write(*,*) 'PREDICTOR DZ'
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
c xmm es la primer matriz de (M,M) en el
!------------------------------------------------------------------------------!
! Codigo del predictor:
!------------------------------------------------------------------------------!
! tdstep predictor es 0.5 tdstep magnum
       tdstep1=timestep*0.5
       write(*,*) 'TDSTEP =', tdstep1
! Paso1: Con las matrices pasadas F1a y F1b extrapolamos a F3----> Extrapolacion
       F3=(7.D0/4.D0)*F1b-(3.D0/4.D0)*F1a
!       F3=1.750D0*F1b-0.750D0*F1a
! Paso2: Usando H3, la matriz densidad rho2 es propagada a rho4----> Prediccion
       rho2t=rho2
       call g2g_timer_start('magnus ins predictor')
       call cumagnusfac(F3,rho2,rho4,M_in,NBCH,tdstep1,factorial)
       call g2g_timer_stop('magnus ins predictor')
! Paso3: Escribimos rho4 en el RMM para poder obtener F5 en el siguiente paso.
!       call g2g_timer_start('cumatmul_predictor')
!       call cumxp(rho4,devPtrX,rho2t,M)
!       call cumpxt(rho2t,devPtrX,rho2t,M)
       call g2g_timer_start('basechange ins predictor')
        rho2t=basechange_cublas(M_in,rho4,devPtrXc,'inv')
       call g2g_timer_stop('basechange ins predictor')
!       call g2g_timer_stop('cumatmul_predictor')
!       call complex_rho_on_to_ao(rho4,devPtrXc,rho2t,M)
!       do j=1,M
!       do k=j,M
!         if(j.eq.k) then
!           RMM(k+(M2-j)*(j-1)/2)=REAL(rho2t(j,k))
!         else
!           RMM(k+(M2-j)*(j-1)/2)=REAL(rho2t(j,k))*2
!         endif
!       enddo
!       enddo
      call sprepack_ctr('L',M,RMM,rho2t(MTB+1:MTB+M,MTB+1:MTB+M))
       DO i=1,M
          DO j=1,M
             if(rho2t(i,j).ne.rho2t(i,j)) stop 'NAN en FBA -predictor'
          ENDDO
       ENDDO
! Paso4: La matriz densidad 4 es usada para calcular F5------> Corrector
      call g2g_timer_start('int3lu + g2g_solve')
      call int3lu(E2)
      call g2g_solve_groups(0,Ex,0)
      call g2g_timer_stop('int3lu + g2g_solve')
      call field_calc(E1, time)
!------------------------------------------------------------------------------!
! Escritura de fock cuadrada
!------------------------------------------------------------------------------!
! Parte Inferior Izquierda (con diagonal)
!       do j=1,M
!         do k=1,j
!           FBA(j,k)=RMM(M5+j+(M2-k)*(k-1)/2-1)
!         enddo
! Parte Superior Derecha (sin diagonal)
!         do k=j+1,M
!           FBA(j,k)=RMM(M5+k+(M2-j)*(j-1)/2-1)
!         enddo
!       enddo

!DFTB: We copy FON inside FBA before this is overwritten to conserve TB terms.
!      This last step is unnecessary if there is not a DFTB calc.

       FBA=FON

       write(*,*) 'FBA escrita'
       call spunpack('L',M,RMM(M5),FBA(MTB+1:MTB+M,MTB+1:MTB+M))
       DO i=1,M
          DO j=1,M
             if(FBA(i,j).ne.FBA(i,j)) stop 'NAN en FBA -predictor'
          ENDDO
       ENDDO
! Ahora tenemos F5 transformada en base de ON y en su forma cuadrada
!       call g2g_timer_start('cumatmul2_predictor')
!       call cumxtf(FBA,devPtrX,FON,M)
!       call cumfx(FON,DevPtrX,FON,M)
!       call g2g_timer_stop('cumatmul2_predictor')
!       call fock_ao_to_on(FBA,devPtrX,FON,M)
       FON=basechange_cublas(M_in, FBA, devPtrX, 'dir')
       DO i=1,M
          DO j=1,M
             if(FON(i,j).ne.FON(i,j)) stop 'NAN en FON -predictor'
          ENDDO
       ENDDO
       DEALLOCATE(rho4,rho2t,F3,FBA)
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine cupredictor_DC(F1a,F1b,FON,rho2,devPtrX,factorial,
     > devPtrXc, timestep, time, M_in, MTB)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!Predictor-Corrector Cheng, V.Vooris.PhysRevB.2006.74.155112
! Esta rutina recibe: F1a,F1b,rho2
! Tira: F5 = F(t+(deltat/2))
       use garcha_mod
       use field_data, only: field
       use field_subs, only: field_calc
       IMPLICIT REAL*8 (a-h,o-z)
       REAL*8,intent(inout) :: F1a(M_in,M_in),F1b(M_in,M_in)
       REAL*8,intent(inout) :: FON(M_in,M_in)
       integer*8,intent(in) :: devPtrX,devPtrXc
      integer, intent(in)    :: M_in
      integer, intent(in)    :: MTB
       REAL*8,allocatable :: F3(:,:),FBA(:,:)
       integer :: i,j,k,kk,stat,M1,M2,MM,M5,M7,M9,MMD,M11,M13,M15,M17
       integer :: M19,M20,M3,M18
       real*8 :: E2, tdstep1, E1
      external CUBLAS_INIT, CUBLAS_SHUTDOWN
      integer CUBLAS_INIT
      REAL*8,intent(in) :: factorial(NBCH), timestep, time
       COMPLEX*8, intent(in) :: rho2(M_in,M_in)
       COMPLEX*8,allocatable :: rho4(:,:),rho2t(:,:)
!-----------------------------------------------------------------------------n
      ALLOCATE(rho4(M_in,M_in),rho2t(M_in,M_in),
     >         F3(M_in,M_in),FBA(M_in,M_in))
      M2=2*M
      MM=M*(M+1)/2
      write(*,*) 'PREDICTOR DC'
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
c xmm es la primer matriz de (M,M) en el
!------------------------------------------------------------------------------!
! Codigo del predictor:
!------------------------------------------------------------------------------!
! tdstep predictor es 0.5 tdstep magnum
       tdstep1=timestep*0.5
! Paso1: Con las matrices pasadas F1a y F1b extrapolamos a F3----> Extrapolacion
       F3=(7.D0/4.D0)*F1b-(3.D0/4.D0)*F1a
!       F3=1.750D0*F1b-0.750D0*F1a
! Paso2: Usando H3, la matriz densidad rho2 es propagada a rho4----> Prediccion
       rho2t=rho2
       call cumagnusfac(F3,rho2,rho4,M_in,NBCH,tdstep1,factorial)
! Paso3: Escribimos rho4 en el RMM para poder obtener F5 en el siguiente paso.
!       call g2g_timer_start('cumatmul_predictor')
!       call cumxp(rho4,devPtrX,rho2t,M)
!       call cumpxt(rho2t,devPtrX,rho2t,M)
!       call g2g_timer_stop('cumatmul_predictor')
!       call complex_rho_on_to_ao(rho4,devPtrXc,rho2t,M)
        rho2t=basechange_cublas(M_in,rho4,devPtrXc,'inv')
!       do j=1,M
!       do k=j,M
!         if(j.eq.k) then
!           RMM(k+(M2-j)*(j-1)/2)=REAL(rho2t(j,k))
!         else
!           RMM(k+(M2-j)*(j-1)/2)=REAL(rho2t(j,k))*2
!         endif
!       enddo
!       enddo
      call sprepack_ctr('L',M,RMM,rho2t(MTB+1:MTB+M,MTB+1:MTB+M))
! Paso4: La matriz densidad 4 es usada para calcular F5------> Corrector
      call int3lu(E2)
      call g2g_solve_groups(0,Ex,0)
      call field_calc(E1, time)
!------------------------------------------------------------------------------!
! Escritura de fock cuadrada
!------------------------------------------------------------------------------!
! Parte Inferior Izquierda (con diagonal)
!       do j=1,M
!         do k=1,j
!           FBA(j,k)=RMM(M5+j+(M2-k)*(k-1)/2-1)
!         enddo
! Parte Superior Derecha (sin diagonal)
!         do k=j+1,M
!           FBA(j,k)=RMM(M5+k+(M2-j)*(j-1)/2-1)
!         enddo
!       enddo

       FBA=FON

       call spunpack('L',M,RMM(M5),FBA(MTB+1:MTB+M,MTB+1:MTB+M))
! Ahora tenemos F5 transformada en base de ON y en su forma cuadrada
!       call g2g_timer_start('cumatmul2_predictor')
!       call cumxtf(FBA,devPtrX,FON,M)
!       call cumfx(FON,DevPtrX,FON,M)
!       call g2g_timer_stop('cumatmul2_predictor')
!       call fock_ao_to_on(FBA,devPtrX,FON,M)
       FON=basechange_cublas(M_in, FBA, devPtrX, 'dir')
       DEALLOCATE(rho4,rho2t,F3,FBA)
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
