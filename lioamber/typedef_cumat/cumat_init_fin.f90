! These subroutines initialise cublas matrices and also destroy them.
! initialise performs both allocation and initialization if needed.
! REAL
subroutine allocate_r(this, size_in, gpu_only, stat)
   implicit none
   integer       , intent(in)              :: size_in
   logical       , intent(in)   , optional :: gpu_only
   class(cumat_r), intent(inout)           :: this
   integer       , intent(inout), optional :: stat

   integer :: my_stat

   my_stat       = 0
   this%mat_size = size_in

   if (this%allocated) return
   if (present(gpu_only)) this%gpu_only = gpu_only

#ifdef CUBLAS
   if (.not. this%gpu_only) allocate(this%matrix(size_in, size_in))

   if (.not. cublas_initialised) call cublas_init()
   my_stat = cublas_alloc(size_in * size_in, REAL_SIZE, &
                          this%cu_pointer)

   if ( my_stat /= 0 ) then
      write(*,'(A)') "  ERROR - CUMATRIX_R%ALLOCATE: CUBLAS memory allocation"&
                     &" failed."
      call cublas_shutdown()
   endif  
#else
   allocate(this%matrix(size_in, size_in))
#endif

   this%allocated = .true.
   if (present(stat)) stat = my_stat
end subroutine allocate_r

subroutine initialise_r(this, size_in, matrix_in, gpu_only)
   implicit none
   integer       , intent(in)           :: size_in
   LIODBLE  , intent(in)           :: matrix_in(size_in, size_in)
   logical       , intent(in), optional :: gpu_only
   class(cumat_r), intent(inout)        :: this

   call this%allocate(size_in, gpu_only)
   call this%set(matrix_in)
end subroutine initialise_r

subroutine destroy_r(this)
   implicit none
   class(cumat_r), intent(inout) :: this

   this%mat_size = 0
   if (allocated(this%matrix)) deallocate(this%matrix)
#ifdef CUBLAS
   call cublas_free(this%cu_pointer)
#endif
end subroutine destroy_r

! COMPLEX
subroutine allocate_x(this, size_in, gpu_only, stat)
   implicit none

   integer       , intent(in)              :: size_in
   logical       , intent(in)   , optional :: gpu_only
   class(cumat_x), intent(inout)           :: this
   integer       , intent(inout), optional :: stat

   integer :: my_stat

   my_stat       = 0
   this%mat_size = size_in

   if (this%allocated) return
   if (present(gpu_only)) this%gpu_only = gpu_only

#ifdef CUBLAS
   if (.not. this%gpu_only) allocate(this%matrix(size_in, size_in))

   if (.not. cublas_initialised) call cublas_init()
   my_stat = cublas_alloc(size_in * size_in, COMPLEX_SIZE, &
                          this%cu_pointer)

   if ( my_stat /= 0 ) then
      write(*,'(A)') "  ERROR - CUMATRIX_X%ALLOCATE: CUBLAS memory allocation"&
                     &" failed."
      call cublas_shutdown()
   endif  
#else
   allocate(this%matrix(size_in, size_in))
#endif

   this%allocated = .true.
   if (present(stat)) stat = my_stat
end subroutine allocate_x

subroutine initialise_x(this, size_in, matrix_in, gpu_only)
   implicit none
   integer       , intent(in)           :: size_in
   TDCOMPLEX     , intent(in)           :: matrix_in(size_in, size_in)
   logical       , intent(in), optional :: gpu_only
   class(cumat_x), intent(inout)        :: this

   call this%allocate(size_in, gpu_only)
   call this%set(matrix_in)

end subroutine initialise_x

subroutine destroy_x(this)
   implicit none
   class(cumat_x), intent(inout) :: this

   this%mat_size = 0
   if (allocated(this%matrix)) deallocate(this%matrix)
#ifdef CUBLAS
   call cublas_free(this%cu_pointer)
#endif
end subroutine destroy_x