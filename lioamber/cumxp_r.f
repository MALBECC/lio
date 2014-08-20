            subroutine cumxp_r(A,devPtrX,C,M)
!!!!!!!!  Hace C=Bt*(A*B) para matrices cuadradas
      implicit none
      integer sizeof_complex
      integer*8 devPtrA
      integer*8,intent(in) :: devPtrX
      integer*8 devPtrScratch1
      integer*8 devPtrScratch2
!      integer*8 devPtrScratch3
!      integer*8 devPtrScratch4
      parameter(sizeof_complex=8)
      integer,intent(in) :: M
      COMPLEX*8 , intent(in) :: A(M,M)
!      REAL*8 , intent(in) :: X(M,M)
      COMPLEX*8, intent(out) :: C(M,M)
      REAL*8 alpha,beta
      REAL*8, dimension (:,:), ALLOCATABLE :: scratch,scratch1,scratch2
      integer i,j,stat
      external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      integer CUBLAS_INIT
      allocate(scratch1(M,M),scratch2(M,M))
      stat=CUBLAS_INIT()
      if (stat.NE.0) then
        write(*,*) "initialization failed -cumxp"
        call CUBLAS_SHUTDOWN
        stop
      endif
      alpha=1.0D0
      beta=0.0D0
!      C=0.0D0
!      write(2344,*), A
      scratch1=REAL(A)
      scratch2=AIMAG(A)
!      write(2345,*), scratch1
!      write(2346,*), scratch2
!--------------------------------------------------------------------!
!      stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPtrX)
      stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPtrScratch1)
      stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPtrScratch2)
!      stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPtrScratch3)
!      stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPtrScratch4)
!--------------------------------------------------------------------!
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_complex,scratch1,M,
     > devPtrScratch1,M)
!      stat = CUBLAS_SET_MATRIX(M,M,sizeof_complex,X,M,devPtrX,M)
!-----------------------REAL-----------------------------------------!
      call CUBLAS_DGEMM ('N','N',M,M,M,alpha,devPtrX
     > ,M ,devPtrScratch1,M, beta, devPtrScratch2,M)
!      call CUBLAS_DGEMM ('N','N',M,M,M,alpha,devPtrScratch
!     > ,M ,devPtrX,M, beta, devPtrScratch1,M)
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_complex,devPtrScratch2,M,
     > scratch1,M)
!----------------------COMPLEX---------------------------------------!
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_complex,scratch2,M,
     > devPtrScratch1,M)
      call CUBLAS_DGEMM ('N','N',M,M,M,alpha,devPtrX
     > ,M ,devPtrScratch1,M, beta, devPtrScratch2,M)
!      call CUBLAS_DGEMM ('N','N',M,M,M,alpha,devPtrScratch
!     > ,M ,devPtrX,M, beta, devPtrScratch2,M)
!--------------------------------------------------------------------!      
!      stat = CUBLAS_GET_MATRIX(M,M,sizeof_complex,devPtrScratch1,M,
!     > scratch1,M)
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_complex,devPtrScratch2,M,
     > scratch2,M)
!--------------------------------------------------------------------!
      C=CMPLX(scratch1,scratch2)
!--------------------------------------------------------------------!
!      call CUBLAS_FREE ( devPtrX )
      call CUBLAS_FREE ( devPtrScratch1 )
      call CUBLAS_FREE ( devPtrScratch2 )
!      call CUBLAS_FREE ( devPtrScratch3 )
!      call CUBLAS_FREE ( devPtrScratch4 )
      call CUBLAS_SHUTDOWN
      DEALLOCATE(scratch1,scratch2)
      return
       end












