            subroutine cumpxt(A,devPtrX,C,M)
!===============================================================================!
!!!!!!!!  Hace C=A*X para matrices cuadradas
!===============================================================================!
      implicit none
      integer sizeof_real
      integer*8 devPtrScratch1
      integer*8 devPtrScratch2
      parameter(sizeof_real=8)
      integer,intent(in) :: M
      integer*8,intent(in) :: devPtrX
      REAL*8 alpha,beta
      REAL*8, dimension (:,:), ALLOCATABLE :: scratch1,scratch2
      integer i,j,stat
      external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      integer CUBLAS_INIT
#ifdef TD_SIMPLE
      COMPLEX*8 , intent(in) :: A(M,M)
      COMPLEX*8, intent(out) :: C(M,M)
#else
      COMPLEX*16 , intent(in) :: A(M,M)
      COMPLEX*16, intent(out) :: C(M,M)
#endif
!-------------------------------------------------------------------------------!
      allocate(scratch1(M,M),scratch2(M,M))
      stat=CUBLAS_INIT()
      if (stat.NE.0) then
        write(*,*) "initialization failed -cumpxt"
        call CUBLAS_SHUTDOWN
        stop
      endif
      alpha=1.0D0
      beta=0.0D0
      scratch1=REAL(A)
      scratch2=AIMAG(A)
!--------------------------------------------------------------------!
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrScratch1)
      if (stat.NE.0) then
        write(*,*) "Allocation failed -cumpxt-1"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrScratch2)
      if (stat.NE.0) then
        write(*,*) "Allocation failed -cumpxt-2"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_real,scratch1,M,
     > devPtrScratch1,M)
      if (stat.NE.0) then
        write(*,*) "Matrix setting failed -cumpxt-1"
        call CUBLAS_SHUTDOWN
        stop
      endif
!-----------------------REAL-----------------------------------------!
      call CUBLAS_DGEMM ('N','T',M,M,M,alpha,devPtrScratch1
     > ,M ,devPtrX,M, beta, devPtrScratch2,M)
      if (stat.NE.0) then
        write(*,*) "DGEMM(1) failed -cumpxt"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_real,devPtrScratch2,M,
     > scratch1,M)
      if (stat.NE.0) then
        write(*,*) "matrix copy failed -cumpxt"
        call CUBLAS_SHUTDOWN
        stop
      endif
!----------------------COMPLEX---------------------------------------!
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_real,scratch2,M,
     > devPtrScratch2,M)
      if (stat.NE.0) then
        write(*,*) "Matrix setting failed -cumpxt"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      call CUBLAS_DGEMM ('N','T',M,M,M,alpha,devPtrScratch2
     > ,M ,devPtrX,M, beta, devPtrScratch1,M)
      if (stat.NE.0) then
        write(*,*) "DGEMM(2) failed -cumpxt"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_real,devPtrScratch1,M,
     > scratch2,M)
      if (stat.NE.0) then
        write(*,*) "matrix copy failed -cumpxt"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      C=CMPLX(scratch1,scratch2)
!--------------------------------------------------------------------!
      call CUBLAS_FREE ( devPtrScratch1 )
      call CUBLAS_FREE ( devPtrScratch2 )
      call CUBLAS_SHUTDOWN
      DEALLOCATE(scratch1,scratch2)
      return
       end

