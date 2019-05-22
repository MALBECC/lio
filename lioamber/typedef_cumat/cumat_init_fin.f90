subroutine initialise_r(this, mat_size_in, matrix_in)
   implicit none
   
   integer       , intent(in)    :: mat_size_in
   real(kind=8)  , intent(in)    :: matrix_in(mat_size_in, mat_size_in)
   class(cumat_r), intent(inout) :: this

   integer :: stat

   this%mat_size = mat_size_in
   allocate(this%matr(mat_size_in, mat_size_in))

   this%matr = matrix_in

#ifdef CUBLAS

   if (.not. cublas_initialised) call cublas_init()

   stat = CUBLAS_ALLOC(this%mat_size * this%mat_size, SIZE_OF_REAL, &
                       this%cu_pointer)
   if ( stat /= 0 ) then
      write(*,'(A)') "  ERROR - CUMATRIX_R%INIT: CUBLAS memory allocation"&
                     &" failed."
      call cublas_shutdown()
      stop
   endif

   stat = CUBLAS_SET_MATRIX(this%mat_size, this%mat_size, SIZE_OF_REAL, &
                            this%matr, this%mat_size, this%cu_pointer,  &
                            this%mat_size)
   if ( stat /= 0 ) then
      write(*,'(A)') "  ERROR - CUMATRIX_R%INIT: CUBLAS set matrix failed."
      call cublas_shutdown()
      stop
   endif

#endif
end subroutine initialise_r

subroutine exterminate_r(this)
   implicit none
   class(cumat_r), intent(inout) :: this

   this%mat_size = 0
   if (allocated(this%matr)) deallocate(this%matr)
   call cublas_free(this%cu_pointer)
end subroutine exterminate_r
