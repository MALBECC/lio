! Grabs an input matrix and multiplies with current cumatrix
! storing the result in mat_result.
subroutine multiply_r(this, output_matrix, input_matrix, mat_size)
#ifdef CUBLAS
   use cublasmath, only: cumxp_r
#endif
   
   implicit none
   class(cumat_r), intent(in)    :: this
   integer       , intent(in)    :: mat_size
   real(kind=8)  , intent(in)    :: input_matrix(:,:)
   real(kind=8)  , intent(inout) :: output_matrix(:,:)

#ifdef CUBLAS
   call cumxp_r(input_matrix, this%cu_pointer, output_matrix, mat_size)
#else
   output_matrix = matmul(this%matrix, input_matrix)
#endif

end subroutine multiply_r

subroutine multiply_x(this, output_matrix, input_matrix, mat_size)
#ifdef CUBLAS
   use cublasmath, only: cumxp
#endif
      
   implicit none
   class(cumat_x), intent(in)    :: this
   integer       , intent(in)    :: mat_size
   TDCOMPLEX     , intent(in)    :: input_matrix(:,:)
   TDCOMPLEX     , intent(inout) :: output_matrix(:,:)

#ifdef CUBLAS
   call cumxp(input_matrix, this%cu_pointer, output_matrix, mat_size)
#else
   output_matrix = matmul(this%matrix, input_matrix)
#endif

end subroutine multiply_x