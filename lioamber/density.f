            subroutine density(M,NCO,X,C)
!=============================================================================!
!!!!!!!!  Construction of the density matrix from the orbital coefficients
!=============================================================================!
      implicit none
      REAL*8, intent(in)   :: X(M,3*M)
      REAL*8, allocatable  :: A(:,:)
      REAL*8,intent(out)   :: C(M,M)
      integer,intent(in)   :: M,NCO
      integer              :: DOSM,i,j
#ifdef CUBLAS
      integer*8 devPtrA,devPtrC
      integer sizeof_real
      parameter(sizeof_real=8)
      integer stat
      external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC,CUBLAS_FREE
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      integer CUBLAS_INIT
#endif
!-----------------------------------------------------------------------------!
      allocate(A(M,NCO))
      DOSM=2*M
      do i=1,M
         do j=1,NCO
            A(i,j)=X(i,DOSM+j)
         enddo
      enddo
      C=0.0D0
#ifdef CUBLAS
! Allocate matrix A on the device
      stat= CUBLAS_ALLOC(M*NCO,sizeof_real,devPtrA)
      if (stat.NE.0) then
        write(*,*) "allocation failed -density"
        call CUBLAS_SHUTDOWN
        stop
      endif
! Allocate matrix C on the device
      stat= CUBLAS_ALLOC(M*M, sizeof_real,devPtrC)
      if (stat.NE.0) then
        write(*,*) "allocation failed -density"
        call CUBLAS_SHUTDOWN
        stop
      endif
! Copy A in the device
      stat = CUBLAS_SET_MATRIX(M,NCO,sizeof_real,A,M,
     > devPtrA,M)
      if (stat.NE.0) then
        write(*,*) "matrix copy to the device failed -density"
        call CUBLAS_SHUTDOWN
        stop 
      endif
! Matrix multiplication
      call CUBLAS_DGEMM ('N','T',M,M,NCO,2.0D0,devPtrA
     > ,M ,devPtrA,M, 0.0D0, devPtrC,M)
! Copy C to in the host
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_real,devPtrC,M,
     > C,M)
! Finalize Variables
      call CUBLAS_FREE(devPtrA)
      call CUBLAS_FREE(devPtrC)
#else
      call DGEMM('N','T',M,M,NCO,2.0D0,A,M,A,M,0.0D0,C,M)
#endif
!-----Only for the close shell case:
      do i=1,M
         do j=1,i-1
            C(i,j)=2.0D0*C(i,j)
         enddo
         do j=i+1,M
            C(i,j)=2.0D0*C(i,j)
         enddo
      enddo
      deallocate(A)
      return
      end





