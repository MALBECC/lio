!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine cupredictor_DZ(F1a,F1b,FON,rho2,devPtrX,factorial,
     >                           devPtrXc, timestep,time,M_in,MTB,dim3)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!Predictor-Corrector Cheng, V.Vooris.PhysRevB.2006.74.155112
! Esta rutina recibe: F1a,F1b,rho2
! Tira: F5 = F(t+(deltat/2))
      use garcha_mod, only: M, Md, open, RMM, nbch, nang, natom, nco,
     > rhoalpha, rhobeta, r, d, Iz, ntatom
      use field_data, only: field
      use field_subs, only: field_calc
      use fockbias_subs , only: fockbias_apply
      use mathsubs
      use general_module
      use faint_cpu, only: int3lu
      implicit none
       integer, intent(in)    :: dim3
       REAL*8,intent(inout)   :: F1a(M_in,M_in,dim3),F1b(M_in,M_in,dim3)
       REAL*8,intent(inout)   :: FON(M_in,M_in,dim3)
       integer*8,intent(in)   :: devPtrX,devPtrXc
       integer, intent(in)    :: M_in
       integer, intent(in)    :: MTB
       integer :: i,j,k,kk,stat, M1,M2,MM,M5,M7,M9,MMD,M11,M13,M15,M17,
     > m18, m19, m20, m3
       REAL*8,allocatable :: F3(:,:,:),FBA(:,:,:)
       real*8 :: E1, E2, tdstep1, ex
      external CUBLAS_INIT, CUBLAS_SHUTDOWN
      integer CUBLAS_INIT
      REAL*8,intent(in) :: factorial(NBCH), timestep, time
       COMPLEX*16, intent(in) :: rho2(M_in,M_in, dim3)
       COMPLEX*16,allocatable :: rho4(:,:,:),rho2t(:,:,:)
!-----------------------------------------------------------------------------n
      ALLOCATE(rho4(M_in,M_in,dim3),rho2t(M_in,M_in,dim3),
     >         F3(M_in,M_in,dim3),FBA(M_in,M_in,dim3))
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
       call cumagnusfac(F3(:,:,1),rho2(:,:,1),rho4(:,:,1),M_in,NBCH,
     >                  tdstep1,factorial)

       if (OPEN) then
         call cumagnusfac(F3(:,:,2),rho2(:,:,2),rho4(:,:,2),M_in,NBCH,
     >                    tdstep1,factorial)
         call g2g_timer_stop('magnus ins predictor')
         call g2g_timer_start('basechange ins predictor')
         rho2t(:,:,1) = basechange_cublas(M_in,rho4(:,:,1),devPtrXc,
     >                  'inv')
         rho2t(:,:,2) = basechange_cublas(M_in,rho4(:,:,2),devPtrXc,
     >                  'inv')
         call g2g_timer_stop('basechange ins predictor')
! Paso3open: Escribimos rho4 en rhoalpha y rhobeta para poder obtener
!            F5 en el siguiente paso.
         call sprepack_ctr('L',M,rhoalpha,
     >                      rho2t(MTB+1:MTB+M,MTB+1:MTB+M,1))
         call sprepack_ctr('L',M,rhobeta,
     >                      rho2t(MTB+1:MTB+M,MTB+1:MTB+M,2))
         RMM(1:MM) = rhoalpha + rhobeta

       else
         call g2g_timer_stop('magnus ins predictor')
         call g2g_timer_start('basechange ins predictor')
         rho2t(:,:,1) = basechange_cublas(M_in,rho4(:,:,1),devPtrXc,
     >                                    'inv')
         call g2g_timer_stop('basechange ins predictor')
! Paso3: Escribimos rho4 en el RMM para poder obtener F5 en el siguiente paso.
         call sprepack_ctr('L',M,RMM,rho2t(MTB+1:MTB+M,MTB+1:MTB+M,1))
       end if

       DO i=1,M
       DO j=1,M
        if(rho2t(i,j,1).ne.rho2t(i,j,1)) stop 'NAN en FBA -predictor'
       ENDDO
       ENDDO

       if(OPEN) then
         DO i=1,M
         DO j=1,M
          if(rho2t(i,j,2).ne.rho2t(i,j,2)) stop 'NAN en FBA -predictor'
         ENDDO
         ENDDO
       end if

! Paso4: La matriz densidad 4 es usada para calcular F5------> Corrector
      call g2g_timer_start('int3lu + g2g_solve')
      call int3lu(E2, RMM(1:MM), RMM(M3:M3+MM), RMM(M5:M5+MM),
     >            RMM(M7:M7+MMd), RMM(M9:M9+MMd), RMM(M11:M11+MM),
     >            open)
      call g2g_solve_groups(0,Ex,0)
      call g2g_timer_stop('int3lu + g2g_solve')
      call field_calc(E1, time, RMM(M3:M3+MM), RMM(M5:M5+MM), r, d,
     > Iz, natom, ntatom, open)

!DFTB: We copy FON inside FBA before this is overwritten to conserve TB terms.
!      This last step is unnecessary if there is not a DFTB calc.

       FBA=FON

       write(*,*) 'FBA escrita'
       call spunpack('L',M,RMM(M5),FBA(MTB+1:MTB+M,MTB+1:MTB+M,1))
!Fockbias:
       call fockbias_apply(time, FBA(MTB+1:MTB+M,MTB+1:MTB+M,1))

       DO i=1,M
       DO j=1,M
          if(FBA(i,j,1).ne.FBA(i,j,1)) stop 'NAN en FBA -predictor'
       ENDDO
       ENDDO

       if (OPEN) then
          call spunpack('L',M,RMM(M3),FBA(MTB+1:MTB+M,MTB+1:MTB+M,2))

!Fockbias:
          call fockbias_apply(time, FBA(MTB+1:MTB+M,MTB+1:MTB+M,2))

          DO i=1,M
          DO j=1,M
             if(FBA(i,j,2).ne.FBA(i,j,2)) stop 'NAN en FBA -predictor'
          ENDDO
          ENDDO
       end if

! Ahora tenemos F5 transformada en base de ON y en su forma cuadrada
       FON(:,:,1)=basechange_cublas(M_in, FBA(:,:,1), devPtrX, 'dir')
       DO i=1,M
       DO j=1,M
          if(FON(i,j,1).ne.FON(i,j,1)) stop 'NAN en FON -predictor'
       ENDDO
       ENDDO

       if(OPEN) then
          FON(:,:,2)=basechange_cublas(M_in, FBA(:,:,2), devPtrX,
     >                                 'dir')
          DO i=1,M
          DO j=1,M
           if(FON(i,j,2).ne.FON(i,j,2)) stop 'NAN en FON -predictor'
          ENDDO
          ENDDO
       end if

       DEALLOCATE(rho4,rho2t,F3,FBA)
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine cupredictor_DC(F1a,F1b,FON,rho2,devPtrX,factorial,
     >                           devPtrXc,timestep,time,M_in,MTB,dim3)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!Predictor-Corrector Cheng, V.Vooris.PhysRevB.2006.74.155112
! Esta rutina recibe: F1a,F1b,rho2
! Tira: F5 = F(t+(deltat/2))
       use garcha_mod, only: M, Md, open, RMM, nbch, nang, natom, nco,
     > rhoalpha, rhobeta, r, d, Iz, ntatom
       use field_data, only: field
       use faint_cpu, only: int3lu
       use field_subs, only: field_calc
       use fockbias_subs , only: fockbias_apply
       implicit none
       integer, intent(in)  :: dim3
       REAL*8,intent(inout) :: F1a(M_in,M_in,dim3),F1b(M_in,M_in,dim3)
       REAL*8,intent(inout) :: FON(M_in,M_in,dim3)
       integer*8,intent(in) :: devPtrX,devPtrXc
      integer, intent(in)    :: M_in
      integer, intent(in)    :: MTB
       REAL*8,allocatable :: F3(:,:,:),FBA(:,:,:)
       integer :: i,j,k,kk,stat,M1,M2,MM,M5,M7,M9,MMD,M11,M13,M15,M17
       integer :: M19,M20,M3,M18
       real*8 :: E2, tdstep1, E1, ex
      external CUBLAS_INIT, CUBLAS_SHUTDOWN
      integer CUBLAS_INIT
      REAL*8,intent(in) :: factorial(NBCH), timestep, time
       COMPLEX*8, intent(in) :: rho2(M_in,M_in,dim3)
       COMPLEX*8,allocatable :: rho4(:,:,:),rho2t(:,:,:)
!-----------------------------------------------------------------------------n
      ALLOCATE(rho4(M_in,M_in,dim3),rho2t(M_in,M_in,dim3),
     >         F3(M_in,M_in,dim3),FBA(M_in,M_in,dim3))
      M2=2*M
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2
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
! Paso2: Usando H3, la matriz densidad rho2 es propagada a rho4----> Prediccion
       rho2t=rho2
       call cumagnusfac(F3(:,:,1),rho2(:,:,1),rho4(:,:,1),M_in,NBCH,
     >                  tdstep1, factorial)
       if (OPEN) then
         call cumagnusfac(F3(:,:,2),rho2(:,:,2),rho4(:,:,2),M_in,NBCH,
     >                    tdstep1,factorial)
         call g2g_timer_stop('magnus ins predictor')
         call g2g_timer_start('basechange ins predictor')
         rho2t(:,:,1) = basechange_cublas(M_in,rho4(:,:,1),devPtrXc,
     >                                    'inv')
         rho2t(:,:,2) = basechange_cublas(M_in,rho4(:,:,2),devPtrXc,
     >                                    'inv')
         call g2g_timer_stop('basechange ins predictor')
! Paso3open: Escribimos rho4 en rhoalpha y rhobeta para poder obtener
!            F5 en el siguiente paso.
         call sprepack_ctr('L',M,rhoalpha,
     >                     rho2t(MTB+1:MTB+M,MTB+1:MTB+M,1))
         call sprepack_ctr('L',M,rhobeta,
     >                     rho2t(MTB+1:MTB+M,MTB+1:MTB+M,2))
         RMM(1:MM) = rhoalpha + rhobeta
       else
         call g2g_timer_stop('magnus ins predictor')
         call g2g_timer_start('basechange ins predictor')
         rho2t(:,:,1) = basechange_cublas(M_in,rho4(:,:,1),devPtrXc,
     >                                    'inv')
         call g2g_timer_stop('basechange ins predictor')
! Paso3: Escribimos rho4 en el RMM para poder obtener F5 en el siguiente paso.
         call sprepack_ctr('L',M,RMM,rho2t(MTB+1:MTB+M,MTB+1:MTB+M,1))
       end if
! Paso4: La matriz densidad 4 es usada para calcular F5------> Corrector
      call int3lu(E2, RMM(1:MM), RMM(M3:M3+MM), RMM(M5:M5+MM),
     >            RMM(M7:M7+MMd), RMM(M9:M9+MMd), RMM(M11:M11+MM),
     >            open)
      call g2g_solve_groups(0,Ex,0)
      call field_calc(E1, time, RMM(M3:M3+MM), RMM(M5:M5+MM), r, d,
     > Iz, natom, ntatom, open)
       FBA=FON

       call spunpack('L',M,RMM(M5),FBA(MTB+1:MTB+M,MTB+1:MTB+M,1))
!Fockbias:
       call fockbias_apply(time, FBA(MTB+1:MTB+M,MTB+1:MTB+M,1))
       FON(:,:,1)=basechange_cublas(M_in, FBA(:,:,1), devPtrX, 'dir')

       if (OPEN) then
          call spunpack('L',M,RMM(M3),FBA(MTB+1:MTB+M,MTB+1:MTB+M,2))
!Fockbias:
          call fockbias_apply(time, FBA(MTB+1:MTB+M,MTB+1:MTB+M,2))
          FON(:,:,2)=basechange_cublas(M_in, FBA(:,:,2), devPtrX,
     >                                 'dir')
       end if

       DEALLOCATE(rho4,rho2t,F3,FBA)
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
