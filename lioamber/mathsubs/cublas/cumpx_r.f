            subroutine cumpx_r(A,devPtrX,C,M)
!!!!!!!!  Hace C=A*X para matrices cuadradas
      implicit none
      integer sizeof_real
      integer*8 devPtrScratch1
      integer*8 devPtrScratch2
      parameter(sizeof_real=8)
      integer,intent(in) :: M
      integer*8,intent(in) :: devPtrX
      REAL*8 , intent(in) :: A(M,M)
      REAL*8, intent(out) :: C(M,M)
      REAL*8 alpha,beta
      REAL*8, dimension (:,:), ALLOCATABLE :: scratch1,scratch2
      integer i,j,stat
      external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      integer CUBLAS_INIT
      allocate(scratch1(M,M),scratch2(M,M))
      alpha=1.0D0
      beta=0.0D0
      scratch1=A
      scratch2=0
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrScratch1)
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrScratch2)
      if (stat.NE.0) then
        write(*,*) "initialization failed -cublas alloc-cumpx_r"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_real,scratch1,M,
     > devPtrScratch1,M)
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_real,scratch2,M,
     > devPtrScratch2,M)
      if (stat.NE.0) then
        write(*,*) "initialization failed -cublas set matrix-cumpx_r"
        call CUBLAS_SHUTDOWN
        stop
      endif
      call CUBLAS_DGEMM ('N','N',M,M,M,alpha,devPtrScratch1
     > ,M ,devPtrX,M, beta, devPtrScratch2,M)
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_real,devPtrScratch2,M,
     > C,M)
      if (stat.NE.0) then
        write(*,*) "initialization failed -cublas get matrix-cumpx_r"
        call CUBLAS_SHUTDOWN
        stop
      endif
      call CUBLAS_FREE ( devPtrScratch1 )
      call CUBLAS_FREE ( devPtrScratch2 )
      DEALLOCATE(scratch1,scratch2)
      return
      end subroutine












