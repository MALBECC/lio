! These subroutines initialise cublas matrices and also destroy them.
subroutine initialise_r(this, mat_size_in, matrix_in)
   implicit none
   
   integer       , intent(in)           :: mat_size_in
   class(cumat_r), intent(inout)        :: this
   real(kind=8)  , intent(in), optional :: matrix_in(mat_size_in, mat_size_in)

   integer :: stat

   this%mat_size = mat_size_in
   allocate(this%matrix(mat_size_in, mat_size_in))

   this%matrix = 0.0D0
   if (present(matrix_in)) this%matrix = matrix_in

#ifdef CUBLAS

   if (.not. cublas_initialised) call cublas_init()

   stat = cublas_alloc(this%mat_size * this%mat_size, SIZE_OF_REAL, &
                       this%cu_pointer)
   if ( stat /= 0 ) then
      write(*,'(A)') "  ERROR - CUMATRIX_R%INIT: CUBLAS memory allocation"&
                     &" failed."
      call cublas_shutdown()
      stop
   endif

   stat = cublas_set_matrix(this%mat_size, this%mat_size, SIZE_OF_REAL,  &
                            this%matrix, this%mat_size, this%cu_pointer, &
                            this%mat_size)
   if ( stat /= 0 ) then
      write(*,'(A)') "  ERROR - CUMATRIX_R%INIT: CUBLAS set matrix failed."
      call cublas_shutdown()
      stop
   endif

#endif
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


subroutine initialise_x(this, mat_size_in, matrix_in)
   implicit none
   integer       , intent(in)           :: mat_size_in
   class(cumat_x), intent(inout)        :: this
   TDCOMPLEX     , intent(in), optional :: matrix_in(mat_size_in, mat_size_in)

   integer :: stat

   this%mat_size = mat_size_in
   allocate(this%matrix(mat_size_in, mat_size_in))

   this%matrix = 0.0D0
   if (present(matrix_in)) this%matrix = matrix_in

#ifdef CUBLAS

   if (.not. cublas_initialised) call cublas_init()

   stat = cublas_alloc(this%mat_size * this%mat_size, SIZE_OF_REAL, &
                       this%cu_pointer)
   if ( stat /= 0 ) then
      write(*,'(A)') "  ERROR - CUMATRIX_R%INIT: CUBLAS memory allocation"&
                     &" failed."
      call cublas_shutdown()
      stop
   endif

   stat = cublas_set_matrix(this%mat_size, this%mat_size, SIZE_OF_REAL,  &
                            this%matrix, this%mat_size, this%cu_pointer, &
                            this%mat_size)
   if ( stat /= 0 ) then
      write(*,'(A)') "  ERROR - CUMATRIX_R%INIT: CUBLAS set matrix failed."
      call cublas_shutdown()
      stop
   endif

#endif
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