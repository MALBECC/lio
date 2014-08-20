      subroutine cuconmut(F,P,C,M)
!!!!!!!!HACE FP-PF
      implicit none
      integer sizeof_complex
      integer*8 devPtrP
      integer*8 devPtrF
      integer*8 devPtrC

      parameter(sizeof_complex=8)
      integer,intent(in) :: M
      COMPLEX*8 , intent(in) :: P(M,M)
      REAL*8 , intent(in) :: F(M,M)
      COMPLEX*8, intent(out) :: C(M,M)
      complex*8 alpha,beta
      COMPLEX*8, dimension (:,:), ALLOCATABLE :: scratch
      integer i,j,stat
      external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      integer CUBLAS_INIT
      stat=CUBLAS_INIT()
      if (stat.NE.0) then
        write(*,*) "initialization failed -cuconmutc"
        call CUBLAS_SHUTDOWN
        stop
      endif
       allocate(scratch(M,M))
       alpha=(1.0000000000,0.00000000000)
       beta=(0.00000000000,0.00000000000)
       C=(0.0D0,0.0D0)
      do i=1,M
      do j=1,M
      scratch(i,j)=CMPLX(F(i,j))
      enddo
      enddo
      stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPtrP)
      stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPtrF)
      stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPtrC)
      if (stat.NE.0) then
        write(*,*) "device memory allocation failed"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_complex,P,M,devPtrP,M)
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_complex,scratch,M,devPtrF,M)
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_complex,C,M,devPtrC,M)
      if (stat.NE.0) then
        call CUBLAS_FREE( devPtrP )
        write(*,*) "data download failed"
        call CUBLAS_SHUTDOWN
          stop
      endif
      call CUBLAS_CGEMM ('N','N',M,M,M,alpha,devPtrF
     > ,M ,devPtrP,M, beta, devPtrC,M)
      beta=(-1.00000000000,0.00000000000)
      call CUBLAS_CGEMM ('CUBLAS_OP_N','CUBLAS_OP_N',M,M,M,alpha,devPtrP
     > ,M ,devPtrF,M, beta, devPtrC,M)
      
      stat = CUBLAS_GET_MATRIX(M, M, sizeof_complex, devPtrC, M, c, M )

      if (stat.NE.0) then
      write(*,*) "data upload failed"
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
