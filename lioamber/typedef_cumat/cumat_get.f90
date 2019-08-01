! These subroutines get matrix values or export matrix values from
! input/output arrays.
! REAL
subroutine get_r(this, matrix_in, stat)
   implicit none
   real(kind=8)  , intent(inout)           :: matrix_in(:,:)
   class(cumat_r), intent(inout)           :: this
   integer       , intent(inout), optional :: stat

   integer :: my_stat
   my_stat = 0

#ifdef CUBLAS
   my_stat = cublas_get_matrix(this%mat_size, this%mat_size, REAL_SIZE,   &
                               this%cu_pointer, this%mat_size, matrix_in, &
                               this%mat_size)
   if ( my_stat /= 0 ) then
      write(*,'(A)') "  ERROR - CUMATRIX_R%GET: CUBLAS get matrix failed."
      call cublas_shutdown()
   endif
   
   ! Updates the cpu version of the matrix. Just because it can.
   if (.not. this%gpu_only ) this%matrix = matrix_in
#else
   matrix_in = this%matrix
#endif

   if (present(stat)) stat = my_stat
end subroutine get_r

subroutine get_rp(this, pointer_in, stat)
   implicit none
   integer(kind=8), intent(in)              :: pointer_in
   class(cumat_r) , intent(inout)           :: this
   integer        , intent(inout), optional :: stat

   integer :: my_stat

   my_stat         = 0
   this%cu_pointer = pointer_in
   
#ifdef CUBLAS
      my_stat = cublas_dcopy(this%mat_size * this%mat_size, this%cu_pointer, &
                             1, pointer_in, 1)
      if ( my_stat /= 0 ) then
         write(*,'(A)') "  ERROR - CUMATRIX_RP%GET: CUBLAS get matrix failed."
         call cublas_shutdown()
      endif
#endif

   if (present(stat)) stat = my_stat
end subroutine get_rp

! COMPLEX
subroutine get_x(this, matrix_in, stat)
   implicit none
   TDCOMPLEX     , intent(inout)           :: matrix_in(:,:)
   class(cumat_x), intent(inout)           :: this
   integer       , intent(inout), optional :: stat

   integer :: my_stat
   my_stat = 0

#ifdef CUBLAS
   my_stat = cublas_get_matrix(this%mat_size, this%mat_size, COMPLEX_SIZE, &
                               this%cu_pointer, this%mat_size, matrix_in,  &
                               this%mat_size)
   if ( my_stat /= 0 ) then
      write(*,'(A)') "  ERROR - CUMATRIX_X%GET: CUBLAS get matrix failed."
      call cublas_shutdown()
   endif

   if (.not. this%gpu_only ) this%matrix = matrix_in
#else
   matrix_in = this%matrix
#endif

   if (present(stat)) stat = my_stat
end subroutine get_x

subroutine get_xp(this, pointer_in, stat)
   implicit none
   integer(kind=8), intent(in)              :: pointer_in
   class(cumat_x) , intent(inout)           :: this
   integer        , intent(inout), optional :: stat

   integer :: my_stat

   my_stat         = 0
   this%cu_pointer = pointer_in
   
#ifdef CUBLAS
   my_stat = cublas_xcopy(this%mat_size * this%mat_size, this%cu_pointer, &
                           1, pointer_in, 1)
   if ( my_stat /= 0 ) then
      write(*,'(A)') "  ERROR - CUMATRIX_XP%GET: CUBLAS get matrix failed."
      call cublas_shutdown()
   endif
#endif

   if (present(stat)) stat = my_stat
end subroutine get_xp