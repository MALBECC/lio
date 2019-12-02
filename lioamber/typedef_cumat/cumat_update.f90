! These subroutines get matrix values or export matrix values from
! input/output arrays.
! REAL
subroutine update_r(this, stat)
   implicit none
   class(cumat_r), intent(inout)           :: this
   integer       , intent(inout), optional :: stat

   integer :: my_stat
   my_stat = 0

#ifdef CUBLAS
   if (this%gpu_only) return
   my_stat = cublas_get_matrix(this%mat_size, this%mat_size, REAL_SIZE,     &
                               this%cu_pointer, this%mat_size, this%matrix, &
                               this%mat_size)
   if ( my_stat /= 0 ) then
      write(*,'(A)') "  ERROR - CUMATRIX_R%UPDATE: CUBLAS get matrix failed."
      call cublas_shutdown()
   endif
#endif

   if (present(stat)) stat = my_stat
end subroutine update_r

! COMPLEX
subroutine update_x(this, stat)
   implicit none
   class(cumat_x), intent(inout)           :: this
   integer       , intent(inout), optional :: stat

   integer :: my_stat
   my_stat = 0

#ifdef CUBLAS
   if (this%gpu_only) return
   my_stat = cublas_get_matrix(this%mat_size, this%mat_size, COMPLEX_SIZE,  &
                               this%cu_pointer, this%mat_size, this%matrix, &
                               this%mat_size)
   if ( my_stat /= 0 ) then
      write(*,'(A)') "  ERROR - CUMATRIX_X%UPDATE: CUBLAS get matrix failed."
      call cublas_shutdown()
   endif
#endif

   if (present(stat)) stat = my_stat
end subroutine update_x