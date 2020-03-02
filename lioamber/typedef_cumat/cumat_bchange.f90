! Changes the base of the input matrix, using the cumat object as the
! basechange matrix. Mode indicates whether the change is direct (dir)
! or inverse (inv).
subroutine change_base_rr(this, input_matrix, mode)
#ifdef CUBLAS
   use cublasmath, only: basechange_cublas
#else
   use mathsubs  , only: basechange_gemm
#endif
   use liosubs_math
   implicit none
   class(cumat_r)  , intent(in)    :: this
   character(len=3), intent(in)    :: mode
   LIODBLE    , intent(inout) :: input_matrix(:,:)

#ifdef CUBLAS
   input_matrix = basechange_cublas(size(input_matrix,1), input_matrix, &
                                    this%cu_pointer, mode)
#else
   input_matrix = basechange_gemm(size(input_matrix,1), input_matrix, &
                                  this%matrix, mode)
#endif

end subroutine change_base_rr

subroutine change_base_xx(this, input_matrix, mode)
#ifdef CUBLAS
   use cublasmath, only: basechange_cublas
#else
   use mathsubs  , only: basechange_gemm
#endif
   implicit none
   class(cumat_x)  , intent(in)    :: this
   character(len=3), intent(in)    :: mode
   TDCOMPLEX       , intent(inout) :: input_matrix(:,:)

#ifdef CUBLAS
   input_matrix = basechange_cublas(size(input_matrix,1), input_matrix, &
                                    this%cu_pointer, mode)
#else
   input_matrix = basechange_gemm(size(input_matrix,1), input_matrix, &
                                  this%matrix, mode)
#endif

end subroutine change_base_xx