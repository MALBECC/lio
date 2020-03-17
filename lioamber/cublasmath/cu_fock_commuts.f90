!-----------------------------------------------------------------------------------------
! Calculates F' and [F',P'] for diis matrices fockm and FP_PFm
!
! Input: F (fock), P (rho), X, Y
!
! F' = X^T * F * X
! P' = Y^T * P * Y
! X = (Y^-1)^T
! => [F',P'] = X^T * F * P * Y - Y^T * P * F * X
! => = A - A^T
! where A = X^T * F * P * Y
!
! Output: A (scratch), A^T (scratch1), F' (fock)
!-----------------------------------------------------------------------------------------
subroutine cu_calc_fock_commuts(fock, rho, devPtrX, devPtrY,commut, M)
   implicit none
   integer        , intent(in)    :: M
   CUDAPTR, intent(in)    :: devPtrX, devPtrY
   LIODBLE   , intent(in)    :: rho(M, M)
   LIODBLE   , intent(out)   :: commut(M, M)
   LIODBLE   , intent(inout) :: fock(M, M)

   LIODBLE, allocatable :: scratch(:,:)
   CUDAPTR :: devPtrFock, devPtrRho, devPtrScr, devPtrScr2
   LIODBLE    :: alpha, beta
   integer         :: i,j,stat, sizeof_real
   parameter(sizeof_real=8)
   external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX, CUBLAS_SHUTDOWN, &
            CUBLAS_ALLOC, CUBLAS_DGEMM, CUBLAS_FREE
   integer  CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX, CUBLAS_INIT, &
            CUBLAS_DGEMM

   alpha = 1.0D0
   beta  = 0.0D0
   allocate(scratch(M, M))  
   ! ALLOCATE STUFF IN THE DEVICE
   stat = CUBLAS_ALLOC(M*M, sizeof_real, devPtrFock)
   stat = CUBLAS_ALLOC(M*M, sizeof_real, devPtrRho) 
   stat = CUBLAS_ALLOC(M*M, sizeof_real, devPtrScr)
   stat = CUBLAS_ALLOC(M*M, sizeof_real, devPtrScr2)

   ! COPY FOCK & RHO TO THE DEVICE
   stat = CUBLAS_SET_MATRIX(M, M, sizeof_real, fock, M, devPtrFock, M)
   stat = CUBLAS_SET_MATRIX(M, M, sizeof_real, rho , M, devPtrRho , M)

   ! X^T * F = scratch^T
   stat = CUBLAS_DGEMM ('T','N', M, M, M, alpha, devPtrX, M, devPtrFock, M, &
                        beta, devPtrScr, M)

   ! do * X for fockm
   stat = CUBLAS_DGEMM ('N','N', M, M, M, alpha, devPtrScr, M, devPtrX, M,  &
                         beta, devPtrScr2, M)

   ! copy back fock matrix to the host
   stat = CUBLAS_GET_MATRIX(M, M, sizeof_real, devPtrScr2, M, fock, M)

   ! * P = scratch1^T
   stat = CUBLAS_DGEMM ('N','N', M, M, M, alpha, devPtrScr, M, devPtrRho, M,&
                        beta, devPtrScr2, M)

   ! * Y = scratch 
   stat = CUBLAS_DGEMM ('N','N', M, M, M, alpha, devPtrScr2, M, devPtrY, M, &
                        beta, devPtrScr, M)

   ! GET THE RESULT BACK TO THE HOST
   stat = CUBLAS_GET_MATRIX(M, M, sizeof_real, devPtrScr, M, scratch, M)

   ! COMMUT = [F',P'] = A - A^T
   do i = 1, M
   do j = 1, M
      commut(i,j)=scratch(i,j)-scratch(j,i)
   enddo
   enddo

   ! Deinitialize stuff
   deallocate(scratch)
   call CUBLAS_FREE(devPtrScr)
   call CUBLAS_FREE(devPtrScr2)
   call CUBLAS_FREE(devPtrRho)
   call CUBLAS_FREE(devPtrFock)
   return
end subroutine
