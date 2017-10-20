            subroutine cumpx(A,devPtrX,C,M)
!========================================================================!
!!!!!!!!  Hace C=A*X para matrices cuadradas
!========================================================================!
      implicit none
      integer*8 sizeof_real
      integer*8 devPtrScratch1
      integer*8 devPtrScratch2
      integer,intent(in) :: M
      integer*8,intent(in) :: devPtrX
      REAL*8 alpha,beta
      REAL*8, dimension (:,:), ALLOCATABLE :: scratch1,scratch2
      integer i,j,stat
      external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC, CUBLAS_DGEMM
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      integer CUBLAS_INIT,CUBLAS_DGEMM
      parameter(sizeof_real=8)
#ifdef TD_SIMPLE
      COMPLEX*8 , intent(in) :: A(M,M)
      COMPLEX*8, intent(out) :: C(M,M)
#else
      COMPLEX*16 , intent(in) :: A(M,M)
      COMPLEX*16, intent(out) :: C(M,M)
#endif
!---------------------------------------------------------------------!
      allocate(scratch1(M,M),scratch2(M,M))
      stat=CUBLAS_INIT()
      if (stat.NE.0) then
        write(*,*) "initialization failed -cumpx"
        call CUBLAS_SHUTDOWN
        stop
      endif
      alpha=1.0D0
      beta=0.0D0
      scratch1=REAL(A)
      write(1000,*) scratch1
      scratch2=AIMAG(A)
      write(1001,*) scratch2
!--------------------------------------------------------------------!
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrScratch1)
      if (stat.NE.0) then
        write(*,*) "Matrix allocation failed -cumpx-1"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrScratch2)
      if (stat.NE.0) then
        write(*,*) "Matrix allocation failed -cumpx-2"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_real,scratch1,M,
     > devPtrScratch1,M)
      if (stat.NE.0) then
        write(*,*) "Matrix setting on device -cumpx-1"
        call CUBLAS_SHUTDOWN
        stop
      endif
!-----------------------REAL-----------------------------------------!
      stat=CUBLAS_DGEMM ('N','N',M,M,M,alpha,devPtrScratch1
     > ,M ,devPtrX,M, beta, devPtrScratch2,M)
      if (stat.NE.0) then
        write(*,*) "GEMM failed -cumpx-1"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_real,devPtrScratch2,M,
     > scratch1,M)
      if (stat.NE.0) then
        write(*,*) "Getting matrix from device failed -cumpx-1"
        call CUBLAS_SHUTDOWN
        stop
      endif
!----------------------COMPLEX---------------------------------------!
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_real,scratch2,M,
     > devPtrScratch2,M)
      if (stat.NE.0) then
        write(*,*) "Matrix setting on device -cumpx-2"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      stat=CUBLAS_DGEMM ('N','N',M,M,M,alpha,devPtrScratch2
     > ,M ,devPtrX,M, beta, devPtrScratch1,M)
      if (stat.NE.0) then
        write(*,*) "GEMM failed -cumpx-2"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!      
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_real,devPtrScratch1,M,
     > scratch2,M)
      if (stat.NE.0) then
        write(*,*) "Getting matrix from device failed -cumpx-2"
        call CUBLAS_SHUTDOWN
        stop
      endif
      C=CMPLX(scratch1,scratch2)
!--------------------------------------------------------------------!
      call CUBLAS_FREE ( devPtrScratch1 )
      call CUBLAS_FREE ( devPtrScratch2 )
      call CUBLAS_SHUTDOWN
      DEALLOCATE(scratch1,scratch2)
      return
       end

