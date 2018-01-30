!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE cumagnusfac_dz(Fock,RhoOld,RhoNew,M,N,dt,factorial)
!------------------------------------------------------------------------------!
! Propagador basado en la formula de Baker-Campbell-Hausdorff de orden N.
! Entrada: Fock(t+(deltat/2)), rho(t)
! Salida:  rho6=rho(t+deltat)
! Ojo:     todo entra y sale en base ortonormal!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       IMPLICIT NONE
       INTEGER,INTENT(IN)                   :: M,N
       REAL*8,INTENT(IN)                    :: Fock(M,M),dt
       INTEGER                              :: ii,stat
       EXTERNAL CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       EXTERNAL CUBLAS_SHUTDOWN, CUBLAS_ALLOC,CUBLAS_CGEMM, CUBLAS_ZCOPY
       EXTERNAL CUBLAS_CAXPY
       INTEGER CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       INTEGER*8 devPOmega
       INTEGER*8 devPPrev
       INTEGER*8 devPNext
       INTEGER*8 devPRho
       INTEGER*8 devPScratch
       INTEGER i,j
       REAL*8 Fact
       REAL*8, INTENT(IN)     :: factorial(N)
       COMPLEX*16 alpha,beta
       COMPLEX*16,ALLOCATABLE,DIMENSION(:,:) :: Omega1
       COMPLEX*16,INTENT(IN)              :: RhoOld(M,M)
       COMPLEX*16,INTENT(OUT)                :: RhoNew(M,M)
       COMPLEX*16,PARAMETER :: icmplx=DCMPLX(0.0D0,1.0D0)
       INTEGER, PARAMETER :: sizeof_complex=16
       INTEGER CUBLAS_ZGEMM,CUBLAS_ZAXPY,CUBLAS_INIT,CUBLAS_ZCOPY
!------------------------------------------------------------------------------!
       ALLOCATE(Omega1(M,M))
! Omega1 (=W) instantiation
       do i=1,M
       do j=1,M
       Omega1(i,j)=(-1)*(icmplx)*(Fock(i,j))*(dt)
       enddo
       enddo
!------------------------------------------------------------------------------!
       stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPOmega)
       stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPPrev)
       stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPNext)
       stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPRho)
!-------------------------------------------------------------------------------!
       stat= CUBLAS_SET_MATRIX(M,M,sizeof_complex,Omega1,M,devPOmega,M)
!=======================================!
      if (stat.NE.0) then
      write(*,*) "data allocation failed -cumagnusfac/1"
      call CUBLAS_SHUTDOWN()
      stop
      endif
!=======================================!
       stat= CUBLAS_SET_MATRIX(M,M,sizeof_complex,RhoOld,M,devPPrev,M)
!=======================================!
      if (stat.NE.0) then
      write(*,*) "data allocation failed -cumagnusfac/2"
      call CUBLAS_SHUTDOWN()
      stop
      endif
!=======================================!
       stat= CUBLAS_SET_MATRIX(M,M,sizeof_complex,RhoOld,M,devPRho,M)
!=======================================!
      if (stat.NE.0) then
      write(*,*) "data allocation failed -cumagnusfac/3"
      call CUBLAS_SHUTDOWN()
      stop
      endif
!=======================================!
! Density matrix propagation
       Rhonew=RhoOld
       alpha=(1.0D0,0.0D0)
       DO ii=1,N
         Fact=factorial(ii)
         alpha=DCMPLX(Fact,0.0D0)
         beta=(0.0D0,0.0D0)
         stat=CUBLAS_ZGEMM('N','N',M,M,M,
     >        alpha,devPOmega,M,devPPrev,M,
     >        beta,devPNext,M)
!=======================================!
      if (stat.NE.0) then
      write(*,*) "CGEM failed -cumagnusfac/1"
      call CUBLAS_FREE ( devPOmega )
      call CUBLAS_FREE ( devPRho )
      call CUBLAS_FREE ( devPNext )
      call CUBLAS_FREE ( devPPrev )
      call CUBLAS_SHUTDOWN()
      stop
      endif
!======================================!
         beta=(-1.0E0,0.0E0)
         stat=CUBLAS_ZGEMM('N','N',M,M,M,
     >        alpha,devPPrev,M,devPOmega,M,
     >        beta,devPNext,M)
!=======================================!
      if (stat.NE.0) then
      write(*,*) "CGEM failed -cumagnusfac/2"
      call CUBLAS_FREE ( devPOmega )
      call CUBLAS_FREE ( devPRho )
      call CUBLAS_FREE ( devPNext )
      call CUBLAS_FREE ( devPPrev )
      endif
!=======================================!
         stat=CUBLAS_ZAXPY(M*M,DCMPLX(1.0D0,0.0D0),devPNext,1,
     >   devPRho,1)
!=======================================!
      if (stat.NE.0) then
      write(*,*) "CAXPY failed -cumagnusfac"
      call CUBLAS_FREE ( devPOmega )
      call CUBLAS_FREE ( devPRho )
      call CUBLAS_FREE ( devPNext )
      call CUBLAS_FREE ( devPPrev )
      call CUBLAS_SHUTDOWN()
      stop
      endif
!=======================================!
      devPScratch=devPPrev
      devPPrev=devPNext
      devPNext=devPScratch
      ENDDO
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_complex,devPRho,M,Rhonew,M)
!=======================================!
      if (stat.NE.0) then
      write(*,*) "data getting failed -cumagnusfac"
      call CUBLAS_SHUTDOWN()
      stop
      endif
!=======================================!
      call CUBLAS_FREE ( devPOmega )
      call CUBLAS_FREE ( devPRho )
      call CUBLAS_FREE ( devPNext )
      call CUBLAS_FREE ( devPPrev )
      DEALLOCATE(Omega1)
      RETURN;END subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE cumagnusfac_dc(Fock,RhoOld,RhoNew,M,N,dt,factorial)
!------------------------------------------------------------------------------!
! Propagador basado en la formula de Baker-Campbell-Hausdorff de orden N.
! Entrada: Fock(t+(deltat/2)), rho(t)
! Salida:  rho6=rho(t+deltat)
! Ojo:     todo entra y sale en base ortonormal!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       IMPLICIT NONE
       INTEGER,INTENT(IN)                   :: M,N
       REAL*8,INTENT(IN)                    :: Fock(M,M),dt
       INTEGER                              :: ii,stat
       EXTERNAL CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       EXTERNAL CUBLAS_SHUTDOWN, CUBLAS_ALLOC,CUBLAS_CGEMM
       EXTERNAL CUBLAS_CAXPY
       INTEGER CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       INTEGER*8 devPOmega
       INTEGER*8 devPPrev
       INTEGER*8 devPNext
       INTEGER*8 devPRho
       INTEGER*8 devPScratch
       INTEGER i,j
       REAL*4 Fact
       REAL*8, INTENT(IN)                   :: factorial(N)
       COMPLEX*8 alpha,beta
       COMPLEX*8,ALLOCATABLE,DIMENSION(:,:) :: Omega1
       COMPLEX*8,INTENT(IN)              :: RhoOld(M,M)
       COMPLEX*8,INTENT(OUT)                :: RhoNew(M,M)
       COMPLEX*8,PARAMETER :: icmplx=CMPLX(0.0D0,1.0D0)
       INTEGER, PARAMETER :: sizeof_complex=8
       INTEGER CUBLAS_CGEMM,CUBLAS_CAXPY,CUBLAS_INIT
!------------------------------------------------------------------------------!
       ALLOCATE(Omega1(M,M))
! Omega1 (=W) instantiation
       do i=1,M
       do j=1,M
       Omega1(i,j)=(icmplx)*(Fock(i,j))*(dt)
       enddo
       enddo
!------------------------------------------------------------------------------!
       stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPOmega)
       stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPPrev)
       stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPNext)
       stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPRho)
!-------------------------------------------------------------------------------!
       stat= CUBLAS_SET_MATRIX(M,M,sizeof_complex,Omega1,M,devPOmega,M)
!=======================================!
      if (stat.NE.0) then
      write(*,*) "data allocation failed -cumagnusfac/1"
      call CUBLAS_SHUTDOWN()
      stop
      endif
!=======================================!
       stat= CUBLAS_SET_MATRIX(M,M,sizeof_complex,RhoOld,M,devPPrev,M)
!=======================================!
      if (stat.NE.0) then
      write(*,*) "data allocation failed -cumagnusfac/2"
      call CUBLAS_SHUTDOWN()
      stop
      endif
!=======================================!
       stat= CUBLAS_SET_MATRIX(M,M,sizeof_complex,RhoOld,M,devPRho,M)
!=======================================!
      if (stat.NE.0) then
      write(*,*) "data allocation failed -cumagnusfac/3"
      call CUBLAS_SHUTDOWN()
      stop
      endif
!=======================================!
! Density matrix propagation
       Rhonew=RhoOld
       alpha=(1.0D0,0.0D0)
       DO ii=1,N
         Fact=factorial(ii)
         alpha=CMPLX(Fact,0.0D0)
         beta=(0.0D0,0.0D0)
         stat=CUBLAS_CGEMM('N','N',M,M,M,
     >        alpha,devPOmega,M,devPPrev,M,
     >        beta,devPNext,M)
!=======================================!
      if (stat.NE.0) then
      write(*,*) "CGEM failed -cumagnusfac/1"
      call CUBLAS_FREE ( devPOmega )
      call CUBLAS_FREE ( devPRho )
      call CUBLAS_FREE ( devPNext )
      call CUBLAS_FREE ( devPPrev )
      call CUBLAS_SHUTDOWN()
      stop
      endif
!======================================!
         beta=(-1.0E0,0.0E0)
         stat=CUBLAS_CGEMM('N','N',M,M,M,
     >        alpha,devPPrev,M,devPOmega,M,
     >        beta,devPNext,M)
!=======================================!
      if (stat.NE.0) then
      write(*,*) "CGEM failed -cumagnusfac/2"
      call CUBLAS_FREE ( devPOmega )
      call CUBLAS_FREE ( devPRho )
      call CUBLAS_FREE ( devPNext )
      call CUBLAS_FREE ( devPPrev )
      call CUBLAS_SHUTDOWN()
      stop
      endif
!=======================================!
         stat=CUBLAS_CAXPY(M*M,cmplx(1.0E0,0.0E0),devPNext,1,
     >   devPRho,1)
!=======================================!
      if (stat.NE.0) then
      write(*,*) "CAXPY failed -cumagnusfac"
      call CUBLAS_FREE ( devPOmega )
      call CUBLAS_FREE ( devPRho )
      call CUBLAS_FREE ( devPNext )
      call CUBLAS_FREE ( devPPrev )
      call CUBLAS_SHUTDOWN()
      stop
      endif
!=======================================!
      devPScratch=devPPrev
      devPPrev=devPNext
      devPNext=devPScratch
      ENDDO
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_complex,devPRho,M,Rhonew,M)
!=======================================!
      if (stat.NE.0) then
      write(*,*) "data getting failed -cumagnusfac"
      call CUBLAS_SHUTDOWN()
      stop
      endif
!=======================================!
      call CUBLAS_FREE ( devPOmega )
      call CUBLAS_FREE ( devPRho )
      call CUBLAS_FREE ( devPNext )
      call CUBLAS_FREE ( devPPrev )
      DEALLOCATE(Omega1)
      RETURN;END subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
