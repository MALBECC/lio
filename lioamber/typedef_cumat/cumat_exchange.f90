! This exchanges the data from another CUMAT with the one from the
! current CUMAT.
! REAL
subroutine exchange_r(this, bmat)
   implicit none
   type(cumat_r) , intent(inout) :: bmat
   class(cumat_r), intent(inout) :: this
   LIODBLE  , allocatable   :: tmp_array(:,:)
#ifdef CUBLAS
   integer(kind=8)               :: tmp_pointer

   tmp_pointer     = bmat%cu_pointer
   bmat%cu_pointer = this%cu_pointer
   this%cu_pointer = tmp_pointer

   if ((.not. this%gpu_only) .and. (.not. bmat%gpu_only)) then
      allocate(tmp_array(this%mat_size, this%mat_size))

      tmp_array   = bmat%matrix
      bmat%matrix = this%matrix
      this%matrix = tmp_array

      deallocate(tmp_array)
   endif

#else
   allocate(tmp_array(this%mat_size, this%mat_size))

   tmp_array   = bmat%matrix
   bmat%matrix = this%matrix
   this%matrix = tmp_array
      
   deallocate(tmp_array)
#endif

end subroutine exchange_r

! COMPLEX
subroutine exchange_x(this, bmat)
   implicit none
   type(cumat_x) , intent(inout) :: bmat
   class(cumat_x), intent(inout) :: this
   TDCOMPLEX     , allocatable   :: tmp_array(:,:)
#ifdef CUBLAS
   integer(kind=8)               :: tmp_pointer

   tmp_pointer     = bmat%cu_pointer
   bmat%cu_pointer = this%cu_pointer
   this%cu_pointer = tmp_pointer

   if ((.not. this%gpu_only) .and. (.not. bmat%gpu_only)) then
      allocate(tmp_array(this%mat_size, this%mat_size))

      tmp_array   = bmat%matrix
      bmat%matrix = this%matrix
      this%matrix = tmp_array

      deallocate(tmp_array)
   endif

#else
   allocate(tmp_array(this%mat_size, this%mat_size))

   tmp_array   = bmat%matrix
   bmat%matrix = this%matrix
   this%matrix = tmp_array
      
   deallocate(tmp_array)
#endif

end subroutine exchange_x
