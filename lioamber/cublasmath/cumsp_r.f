            subroutine cumsp_r(A,devPtrS,C,M)
!!!!!!!!  Hace C=S*A para matrices cuadradas
      implicit none
      integer sizeof_real
      integer*8 devPtrA
      integer*8 devPtrC
      integer*8,intent(in) :: devPtrS
      parameter(sizeof_real=8)
      integer,intent(in) :: M
      REAL*8 , intent(in) :: A(M,M)
      REAL*8, intent(out) :: C(M,M)
      REAL*8 alpha,beta
      integer i,j,stat
      external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC,CUBLAS_DGEMM
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      integer CUBLAS_INIT,CUBLAS_DGEMM
      stat=CUBLAS_INIT()
      if (stat.NE.0) then
        write(*,*) "initialization failed -cumsp_r"
        call CUBLAS_SHUTDOWN
        stop
      endif
      alpha=1.0D0
      beta=0.0D0
!--------------------------------------------------------------------!
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrA)
      if (stat.NE.0) then
        write(*,*) "allocation failed -cumsp_r - 1 -"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrC)
      if (stat.NE.0) then
        write(*,*) "allocation failed -cumsp_r - 3 -"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_real,A,M,devPtrA,M)
      if (stat.NE.0) then
        write(*,*) "matrix setting failed -cumsp_r"
        call CUBLAS_SHUTDOWN
        stop
      endif
!-----------------------GEMM-----------------------------------------!
      stat=CUBLAS_DGEMM ('N','N',M,M,M,alpha,devPtrA,M,
     > devPtrS,M,beta, devPtrC,M)
      if (stat.NE.0) then
        write(*,*) "CUBLAS DGEMM failed -cumsp_r"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_real,devPtrC,M,C,M)
      if (stat.NE.0) then
        write(*,*) "matrix getting failed -cumsp_r"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      call CUBLAS_FREE ( devPtrA )
      call CUBLAS_FREE ( devPtrC )
      call CUBLAS_SHUTDOWN
      return
      end subroutine
