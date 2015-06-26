            subroutine cumxtf(A,devPtrX,C,M)
!!!!!!!!  Hace C=Bt*(A*B) para matrices cuadradas
      implicit none
      integer sizeof_real
      integer*8 devPtrA
      integer*8,intent(in) :: devPtrX
      integer*8 devPtrScratch1
      integer*8 devPtrScratch2
      parameter(sizeof_real=8)
      integer,intent(in) :: M
      REAL*8 , intent(in) :: A(M,M)
      REAL*8, intent(out) :: C(M,M)
      REAL*8 alpha,beta
      integer i,j,stat
      external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      integer CUBLAS_INIT
      alpha=1.0D0
      beta=0.0D0
!--------------------------------------------------------------------!
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrScratch1)
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrScratch2)
!--------------------------------------------------------------------!
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_real,A,M,devPtrScratch1,M)
!-----------------------REAL-----------------------------------------!
      call CUBLAS_DGEMM ('T','N',M,M,M,alpha,devPtrX,M,devPtrScratch1,M,
     > beta, devPtrScratch2,M)
!--------------------------------------------------------------------!
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_real,devPtrScratch2,M,C,M)
!--------------------------------------------------------------------!
      call CUBLAS_FREE ( devPtrScratch1 )
      call CUBLAS_FREE ( devPtrScratch2 )
      return
      end subroutine












