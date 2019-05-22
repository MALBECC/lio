subroutine multiply_r(this, mat_result, mat_inp, mat_size)
   use cublasmath, only: cumxp_r
   
   implicit none
   class(cumat_r), intent(in)    :: this
   integer       , intent(in)    :: mat_size
   real(kind=8)  , intent(in)    :: mat_inp(:,:)
   real(kind=8)  , intent(inout) :: mat_result(:,:)

#ifdef CUBLAS
   call cumxp_r(mat_inp, this%cu_pointer, mat_result, mat_size)
#else
   mat_result = matmul(this%matrix, mat_inp)
#endif

end subroutine multiply_r