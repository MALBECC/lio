      subroutine cuconmut_r(F,P,C,M)
!!!!!!!!HACE FP-PF
      implicit none
      integer sizeof_real
      integer*8 devPtrP
      integer*8 devPtrF
      integer*8 devPtrC
      parameter(sizeof_real=8)
      integer,intent(in) :: M
      REAL*8 , intent(in) :: P(M,M)
      REAL*8 , intent(in) :: F(M,M)
      REAL*8, intent(out) :: C(M,M)
      REAL*8 alpha,beta
      REAL*8, dimension (:,:), ALLOCATABLE :: scratch
      integer i,j,stat
      external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      integer CUBLAS_INIT
      stat=CUBLAS_INIT()
      if (stat.NE.0) then
        write(*,*) "initialization failed -cuconmutc_r"
        call CUBLAS_SHUTDOWN
        stop
      endif
       allocate(scratch(M,M))
       alpha=1.0000000000
       beta=0.00000000000
       C=0.0D0
      do i=1,M
      do j=1,M
      scratch(i,j)=F(i,j)
      enddo
      enddo
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrP)
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrF)
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrC)
      if (stat.NE.0) then
        write(*,*) "device memory allocation failed -cuconmutc_r"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_real,P,M,devPtrP,M)
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_real,scratch,M,devPtrF,M)
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_real,C,M,devPtrC,M)
      if (stat.NE.0) then
        call CUBLAS_FREE( devPtrP )
        write(*,*) "data download failed -cuconmutc_r"
        call CUBLAS_SHUTDOWN
        stop
      endif
      call CUBLAS_DGEMM ('N','N',M,M,M,alpha,devPtrF
     > ,M ,devPtrP,M, beta, devPtrC,M)
      beta=(-1.00000000000,0.00000000000)
      call CUBLAS_DGEMM ('CUBLAS_OP_N','CUBLAS_OP_N',M,M,M,alpha,devPtrP
     > ,M ,devPtrF,M, beta, devPtrC,M)      
      stat = CUBLAS_GET_MATRIX(M, M, sizeof_real, devPtrC, M, c, M )
      if (stat.NE.0) then
      write(*,*) "data upload failed -cuconmutc_r"
      call CUBLAS_FREE ( devPtrP )
      call CUBLAS_FREE ( devPtrF )
      call CUBLAS_FREE ( devPtrC )
      call CUBLAS_SHUTDOWN
      stop
      endif
      call CUBLAS_FREE ( devPtrP )
      call CUBLAS_FREE ( devPtrF )
      call CUBLAS_FREE ( devPtrC )
      call CUBLAS_SHUTDOWN
      DEALLOCATE(scratch)
      end
