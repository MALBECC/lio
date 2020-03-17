! Includes miscellaneous routines to get cumat properties of a given
! matrix (size, allocation, pointers and the like).
! REAL and COMPLEX
function gpu_pointer(this) result(ptr_out)
   implicit none
   class(cumat_basic), intent(in) :: this
   CUDAPTR :: ptr_out

   ptr_out = this%cu_pointer
end function gpu_pointer

function get_size(this) result(size_out)
   implicit none
   class(cumat_basic), intent(in) :: this
   integer :: size_out

   size_out = this%mat_size
end function get_size

function is_allocated(this) result(is_alloc)
   implicit none
   class(cumat_basic), intent(in) :: this
   logical :: is_alloc

   is_alloc = this%allocated
end function is_allocated

function is_gpu_only(this) result(is_gpu)
   implicit none
   class(cumat_basic), intent(in)  :: this
   logical :: is_gpu

   is_gpu = this%gpu_only
end function is_gpu_only