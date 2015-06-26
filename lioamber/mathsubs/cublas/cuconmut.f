!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!     CONMUTATOR SUBROUTINE - cublas version -
!     calculates the commutator [F,P] where F and P are fock and density matrix respectively.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      subroutine cuconmut_dz(F,P,C,M)
!-------------------------------------------------------------------!
      implicit none
      integer sizeof_complex
      integer*8 devPtrP
      integer*8 devPtrF
      integer*8 devPtrC
      integer,intent(in) :: M
      REAL*8 , intent(in) :: F(M,M)
      integer i,j,stat
      external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      integer CUBLAS_INIT
      parameter(sizeof_complex=16)
      COMPLEX*16 , intent(in) :: P(M,M)
      COMPLEX*16, intent(out) :: C(M,M)
      complex*16 :: alpha,beta
      COMPLEX*16, dimension (:,:), ALLOCATABLE :: scratch
!-------------------------------------------------------------------!
       allocate(scratch(M,M))
       alpha=(1,0)
       beta=(0,0)
       C=(0,0)
      do i=1,M
      do j=1,M
      scratch(i,j)=CMPLX(F(i,j),0)
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
!
      call CUBLAS_ZGEMM ('N','N',M,M,M,alpha,devPtrF
     > ,M ,devPtrP,M, beta, devPtrC,M)
      beta=(-1,0)
      call CUBLAS_ZGEMM ('N','N',M,M,M,alpha,
     > devPtrP,M ,devPtrF,M, beta, devPtrC,M)
!
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
      DEALLOCATE(scratch)
      RETURN;END SUBROUTINE
!================================================================================!
      subroutine cuconmut_dc(F,P,C,M)
!--------------------------------------------------------------------------------!
      implicit none
      integer sizeof_complex
      integer*8 devPtrP
      integer*8 devPtrF
      integer*8 devPtrC
      integer,intent(in) :: M
      REAL*8 , intent(in) :: F(M,M)
      integer i,j,stat
      external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      integer CUBLAS_INIT
      parameter(sizeof_complex=8)
      COMPLEX*8 , intent(in) :: P(M,M)
      COMPLEX*8, intent(out) :: C(M,M)
      complex*8 :: alpha,beta
      COMPLEX*8, dimension (:,:), ALLOCATABLE :: scratch
!--------------------------------------------------------------------------------!
       allocate(scratch(M,M))
       alpha=(1,0)
       beta=(0,0)
       C=(0,0)
      do i=1,M
      do j=1,M
      scratch(i,j)=CMPLX(F(i,j),0)
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
!
      call CUBLAS_CGEMM ('N','N',M,M,M,alpha,devPtrF
     > ,M ,devPtrP,M, beta, devPtrC,M)
      beta=(-1,0)
      call CUBLAS_CGEMM ('N','N',M,M,M,alpha,devPtrP
     > ,M ,devPtrF,M, beta, devPtrC,M)
!
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
      DEALLOCATE(scratch)
      RETURN;END SUBROUTINE
!==========================================================================!
      
