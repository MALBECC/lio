!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!   Multiplicates two real matrices F and X
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine cumfx(A, devPtrX,C, M)
   implicit none
   integer        , intent(in)  :: M
   integer(kind=8), intent(in)  :: devPtrX
   LIODBLE   , intent(in)  :: A(M, M)
   LIODBLE   , intent(out) :: C(M, M)

   integer(kind=8) :: devPtrScratch1, devPtrScratch2
   LIODBLE    :: alpha, beta
   integer         :: stat, sizeof_real
   parameter(sizeof_real=8)

   external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX, CUBLAS_SHUTDOWN,&
            CUBLAS_ALLOC, CUBLAS_DGEMM
   integer  CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX, CUBLAS_INIT, &
            CUBLAS_DGEMM

   alpha = 1.0D0
   beta  = 0.0D0
   stat = CUBLAS_ALLOC(M*M, sizeof_real, devPtrScratch1)
   stat = CUBLAS_ALLOC(M*M, sizeof_real, devPtrScratch2)

   stat = CUBLAS_SET_MATRIX(M, M, sizeof_real,A, M, devPtrScratch1, M)
   stat = CUBLAS_DGEMM ('N','N', M, M, M, alpha, devPtrScratch1, M, devPtrX, M,&
                        beta, devPtrScratch2, M)

   stat = CUBLAS_GET_MATRIX(M, M, sizeof_real, devPtrScratch2, M, C, M)

   call CUBLAS_FREE ( devPtrScratch1 )
   call CUBLAS_FREE ( devPtrScratch2 )
   return
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!











