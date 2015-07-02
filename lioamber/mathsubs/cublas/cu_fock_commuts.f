!-----------------------------------------------------------------------------------------
! Calculates F' and [F',P'] for diis matrices fockm and FP_PFm
!
! Input: F (fock), P (rho), X, Y
!
! F' = X^T * F * X
! P' = Y^T * P * Y
! X = (Y^-1)^T
! => [F',P'] = X^T * F * P * Y - Y^T * P * F * X
! =>         = A - A^T
! where A = X^T * F * P * Y
!
! Output: A (scratch), A^T (scratch1), F' (fock)
!-----------------------------------------------------------------------------------------
      subroutine cu_calc_fock_commuts(fock,rho,devPtrX,devPtrY,commut,M)
          
          implicit none
          integer, intent(in)    :: M
          REAL*8,  intent(in)    :: rho(M,M)
          REAL*8,  intent(inout) :: fock(M,M)
          REAL*8,  intent(out)   :: commut(M,M)
          integer sizeof_real
          integer*8 devPtrFock,devPtrRho
          integer*8,intent(in) :: devPtrX,devPtrY
          integer*8 devPtrScratch,devPtrScratch2
          parameter(sizeof_real=8)
          REAL*8 alpha,beta
          REAL*8, allocatable :: scratch(:,:)
          integer i,j,stat
          external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
          external CUBLAS_SHUTDOWN, CUBLAS_ALLOC,CUBLAS_DGEMM
          external CUBLAS_FREE
          integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
          integer CUBLAS_INIT,CUBLAS_DGEMM
!---------------------------------------------------------------------
!         TESTING
!---------------------------------------------------------------------
!          REAL*8, intent(in) :: x(M,M),y(M,M)
!          REAL*8,allocatable :: scratch1(:,:)
!          allocate(scratch1(M,M)) 
!---------------------------------------------------------------------
          alpha=1.0D0
          beta=0.0D0
          allocate(scratch(M,M))  
!---------------------------------------------------------------------
! ALLOCATE STUFF IN THE DEVICE
!---------------------------------------------------------------------
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrFock)
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrRho) 
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrScratch)
      stat= CUBLAS_ALLOC(M*M, sizeof_real, devPtrScratch2)
!---------------------------------------------------------------------
! COPY FOCK & RHO TO THE DEVICE
!---------------------------------------------------------------------
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_real,fock,M,devPtrFock,M)
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_real,rho,M,devPtrRho,M)
!---------------------------------------------------------------------
! X^T * F = scratch^T
!---------------------------------------------------------------------
         stat = CUBLAS_DGEMM ('T','N',M,M,M,alpha,devPtrX,M,
     > devPtrFock,M, beta, devPtrScratch,M)
!---------------------------------------------------------------------
! do * X for fockm
!---------------------------------------------------------------------
         stat = CUBLAS_DGEMM ('N','N',M,M,M,alpha,devPtrScratch,M,
     > devPtrX,M,beta, devPtrScratch2,M)
!---------------------------------------------------------------------
! copy back fock matrix to the host
!---------------------------------------------------------------------
         stat = CUBLAS_GET_MATRIX(M,M,sizeof_real,devPtrScratch2,M,
     > fock,M)        
!---------------------------------------------------------------------
! * P = scratch1^T
!---------------------------------------------------------------------
          stat = CUBLAS_DGEMM ('N','N',M,M,M,alpha,devPtrScratch,M,
     > devPtrRho,M,beta, devPtrScratch2,M)
!---------------------------------------------------------------------
! * Y = scratch 
!---------------------------------------------------------------------
          stat = CUBLAS_DGEMM ('N','N',M,M,M,alpha,devPtrScratch2,M,
     > devPtrY,M,beta, devPtrScratch,M)
!---------------------------------------------------------------------
!  GET THE RESULT BACK TO THE HOST
!---------------------------------------------------------------------
          stat = CUBLAS_GET_MATRIX(M,M,sizeof_real,devPtrScratch,M,
     > scratch,M)
!---------------------------------------------------------------------
!   COMMUT = [F',P'] = A - A^T
!---------------------------------------------------------------------
      DO i=1,M
         DO j=1,M
            commut(i,j)=scratch(i,j)-scratch(j,i)
         ENDDO
      ENDDO
!---------------------------------------------------------------------
! Deinitialize stuff
!---------------------------------------------------------------------
      DEALLOCATE(scratch)
      call CUBLAS_FREE(devPtrScratch)
      call CUBLAS_FREE(devPtrScratch2)
      call CUBLAS_FREE(devPtrRho)
      call CUBLAS_FREE(devPtrFock)
      return
      end subroutine
