            subroutine cumxtp(A,devPtrX,C,M)
!===============================================================================!
!!!!!!!!  Hace C=Bt*(A*B) para matrices cuadradas
!===============================================================================!
      implicit none
      integer sizeof_real
      integer*8 devPtrA
      integer*8,intent(in) :: devPtrX
      integer*8 devPtrScratch1
      integer*8 devPtrScratch2
      parameter(sizeof_real=8)
      integer,intent(in) :: M
      REAL*8 alpha,beta
      REAL*8, dimension (:,:), ALLOCATABLE :: scratch,scratch1,scratch2
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
!---------------------------------------------------------------------!
      allocate(scratch1(M,M),scratch2(M,M))
      stat=CUBLAS_INIT()
      if (stat.NE.0) then
        write(*,*) "initialization failed -cumxtp"
        call CUBLAS_SHUTDOWN
        stop
      endif
      alpha=1.0D0
      beta=0.0D0
      scratch1=REAL(A)
      scratch2=AIMAG(A)
!--------------------------------------------------------------------!
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrScratch1)
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrScratch2)
!--------------------------------------------------------------------!
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_real,scratch1,M,
     > devPtrScratch1,M)
!------------------------REAL----------------------------------------!
      call CUBLAS_DGEMM ('T','N',M,M,M,alpha,devPtrX
     > ,M ,devPtrScratch1,M, beta, devPtrScratch2,M)
!--------------------------------------------------------------------!
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_real,devPtrScratch2,M,
     > scratch1,M)
!----------------------COMPLEX---------------------------------------!
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_real,scratch2,M,
     > devPtrScratch1,M)
!--------------------------------------------------------------------!
      call CUBLAS_DGEMM ('T','N',M,M,M,alpha,devPtrX
     > ,M ,devPtrScratch1,M, beta, devPtrScratch2,M)
!--------------------------------------------------------------------!      
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_real,devPtrScratch2,M,
     > scratch2,M)
!--------------------------------------------------------------------!
      C=CMPLX(scratch1,scratch2)
!--------------------------------------------------------------------!
      call CUBLAS_FREE ( devPtrScratch1 )
      call CUBLAS_FREE ( devPtrScratch2 )
      call CUBLAS_SHUTDOWN
      DEALLOCATE(scratch1,scratch2)
      return
      end












