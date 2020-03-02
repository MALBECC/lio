! Grabs an input matrix (CPU) and multiplies with current cumat
! storing the result in mat_result.
subroutine multiply_r(this, output_matrix, input_matrix)
#ifdef CUBLAS
   use cublasmath, only: cumxp_r
#endif
   
   implicit none
   class(cumat_r), intent(in)    :: this
   LIODBLE  , intent(in)    :: input_matrix(:,:)
   LIODBLE  , intent(inout) :: output_matrix(:,:)

   integer :: mat_size1, mat_size2
   
   mat_size1 = size(input_matrix,1)
   mat_size2 = size(input_matrix,2)
   
#ifdef CUBLAS
   call cumxp_r(input_matrix, this%cu_pointer, output_matrix, mat_size1)
#else
   output_matrix = matmul(this%matrix, input_matrix)
#endif

end subroutine multiply_r

subroutine multiply_x(this, output_matrix, input_matrix)
#ifdef CUBLAS
   use cublasmath, only: cumxp
#endif
      
   implicit none
   class(cumat_x), intent(in)    :: this
   TDCOMPLEX     , intent(in)    :: input_matrix(:,:)
   TDCOMPLEX     , intent(inout) :: output_matrix(:,:)

   integer :: mat_size1, mat_size2
   
   mat_size1 = size(input_matrix,1)
   mat_size2 = size(input_matrix,2)

#ifdef CUBLAS
   call cumxp(input_matrix, this%cu_pointer, output_matrix, mat_size1)
#else
   output_matrix = matmul(this%matrix, input_matrix)
#endif

end subroutine multiply_x