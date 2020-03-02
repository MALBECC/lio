! These subroutines set matrix values or export matrix values from
! input/output arrays.
! REAL
subroutine set_r(this, matrix_in, stat)
   implicit none
   LIODBLE  , intent(in)              :: matrix_in(:,:)
   class(cumat_r), intent(inout)           :: this
   integer       , intent(inout), optional :: stat

   integer :: my_stat
   my_stat = 0

   if (.not. this%allocated) call this%allocate(size(matrix_in,1))

#ifdef CUBLAS
   my_stat = cublas_set_matrix(this%mat_size, this%mat_size, REAL_SIZE,   &
                               matrix_in, this%mat_size, this%cu_pointer, &
                               this%mat_size)
   if ( my_stat /= 0 ) then
      write(*,'(A)') "  ERROR - CUMATRIX_R%SET: CUBLAS set matrix failed."
      call cublas_shutdown()
   endif

   if (.not. this%gpu_only ) this%matrix = matrix_in
#else
   this%matrix = matrix_in
#endif

   if (present(stat)) stat = my_stat
end subroutine set_r

subroutine set_rp(this, pointer_in, size_in, gpu_only_in, stat)
   implicit none
   integer(kind=8), intent(in)              :: pointer_in
   integer        , intent(in)              :: size_in
   logical        , intent(in)   , optional :: gpu_only_in
   class(cumat_r) , intent(inout)           :: this
   integer        , intent(inout), optional :: stat

   integer :: my_stat

   my_stat         = 0
   this%cu_pointer = pointer_in
   this%mat_size   = size_in
   this%gpu_only   = .true.
   if (present(gpu_only_in)) this%gpu_only = gpu_only_in

#ifdef CUBLAS
   if (.not. this%gpu_only) then
      if (.not. allocated(this%matrix)) allocate(this%matrix(size_in,size_in))
      my_stat = cublas_get_matrix(size_in, size_in, REAL_SIZE,           &
                                  this%cu_pointer, size_in, this%matrix, &
                                  size_in)
      if ( my_stat /= 0 ) then
         write(*,'(A)') "  ERROR - CUMATRIX_RP%SET: CUBLAS get matrix failed."
         call cublas_shutdown()
      endif
   endif
#endif

   if (present(stat)) stat = my_stat
end subroutine set_rp

! COMPLEX
subroutine set_x(this, matrix_in, stat)
   implicit none
   TDCOMPLEX     , intent(in)              :: matrix_in(:,:)
   class(cumat_x), intent(inout)           :: this
   integer       , intent(inout), optional :: stat

   integer :: my_stat
   my_stat = 0

   if (.not. this%allocated) call this%allocate(size(matrix_in,1))

#ifdef CUBLAS
   my_stat = cublas_set_matrix(this%mat_size, this%mat_size, COMPLEX_SIZE, &
                               matrix_in, this%mat_size, this%cu_pointer,  &
                               this%mat_size)
   if ( my_stat /= 0 ) then
      write(*,'(A)') "  ERROR - CUMATRIX_X%SET: CUBLAS set matrix failed."
      call cublas_shutdown()
   endif

   if (.not. this%gpu_only ) this%matrix = matrix_in
#else
   this%matrix = matrix_in
#endif

   if (present(stat)) stat = my_stat
end subroutine set_x

subroutine set_xp(this, pointer_in, size_in, gpu_only_in, stat)
   implicit none
   integer(kind=8), intent(in)              :: pointer_in
   integer        , intent(in)              :: size_in
   logical        , intent(in)   , optional :: gpu_only_in
   class(cumat_x) , intent(inout)           :: this
   integer        , intent(inout), optional :: stat

   integer :: my_stat

   my_stat         = 0
   this%cu_pointer = pointer_in
   this%mat_size   = size_in
   this%gpu_only   = .true.
   if (present(gpu_only_in)) this%gpu_only = gpu_only_in

#ifdef CUBLAS
   if (.not. this%gpu_only) then
      if (.not. allocated(this%matrix)) allocate(this%matrix(size_in,size_in))
      my_stat = cublas_get_matrix(size_in, size_in, COMPLEX_SIZE,        &
                                  this%cu_pointer, size_in, this%matrix, &
                                  size_in)
      if ( my_stat /= 0 ) then
         write(*,'(A)') "  ERROR - CUMATRIX_XP%SET: CUBLAS get matrix failed."
         call cublas_shutdown()
      endif
   endif
#endif

   if (present(stat)) stat = my_stat
end subroutine set_xp