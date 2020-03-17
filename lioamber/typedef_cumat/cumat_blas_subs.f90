! This file contains masks for BLAS and CUBLAS subroutines to be performed
! with cumats. If possible, all operations are performed in GPU. The calling
! object is the one to store the result.
!
! * add_mat performs A = a*B + A; B being an external CUMAT (see BLAS AXPY)
! * mat_mul performs A = a*B*C + b*A; B,C being external CUMATs (see BLAS GEMM)

! REAL
subroutine add_mat_r(this, bmat, alfa, stat)
   implicit none
   LIODBLE  , intent(in)              :: alfa
   type(cumat_r) , intent(in)              :: bmat
   class(cumat_r), intent(inout)           :: this
   integer       , intent(inout), optional :: stat

   integer :: my_stat = 0

#ifdef CUBLAS
   my_stat = cublas_daxpy(this%mat_size * this%mat_size, alfa, &
                          bmat%cu_pointer, 1, this%cu_pointer, 1)
#else
   call daxpy(this%mat_size * this%mat_size, alfa, bmat%matrix, 1, &
              this%matrix, 1)
#endif

   if (present(stat)) stat = my_stat
end subroutine add_mat_r

subroutine mat_mul_r(this, bmat, cmat, alfa, beta, stat)
   implicit none
   LIODBLE  , intent(in)              :: alfa, beta
   type(cumat_r) , intent(in)              :: bmat, cmat
   class(cumat_r), intent(inout)           :: this
   integer       , intent(inout), optional :: stat

   integer :: my_stat = 0

#ifdef CUBLAS
   my_stat = cublas_dgemm('N', 'N', this%mat_size, this%mat_size,        &
                          this%mat_size, alfa, bmat%cu_pointer,          &
                          bmat%mat_size, cmat%cu_pointer, cmat%mat_size, &
                          beta, this%cu_pointer, this%mat_size)
#else
   call dgemm('N', 'N', this%mat_size, this%mat_size, this%mat_size, alfa, &
              bmat%matrix, bmat%mat_size, cmat%matrix, cmat%mat_size, beta,&
              this%matrix, this%mat_size)
#endif

   if (present(stat)) stat = my_stat
end subroutine mat_mul_r

! COMPLEX
subroutine add_mat_x(this, bmat, alfa, stat)
   implicit none
   TDCOMPLEX     , intent(in)              :: alfa
   type(cumat_x) , intent(in)              :: bmat
   class(cumat_x), intent(inout)           :: this
   integer       , intent(inout), optional :: stat

   integer :: my_stat = 0

#ifdef CUBLAS
   my_stat = cublas_xaxpy(this%mat_size * this%mat_size, alfa, &
                          bmat%cu_pointer, 1, this%cu_pointer, 1)
#else
   call xaxpy(this%mat_size * this%mat_size, alfa, bmat%matrix, 1, &
              this%matrix, 1)
#endif

   if (present(stat)) stat = my_stat
end subroutine add_mat_x

subroutine mat_mul_x(this, bmat, cmat, alfa, beta, stat)
   implicit none
   TDCOMPLEX     , intent(in)              :: alfa, beta
   type(cumat_x) , intent(in)              :: bmat, cmat
   class(cumat_x), intent(inout)           :: this
   integer       , intent(inout), optional :: stat

   integer :: my_stat = 0

#ifdef CUBLAS
   my_stat = cublas_xgemm('N', 'N', this%mat_size, this%mat_size,        &
                          this%mat_size, alfa, bmat%cu_pointer,          &
                          bmat%mat_size, cmat%cu_pointer, cmat%mat_size, &
                          beta, this%cu_pointer, this%mat_size)
#else
   call xgemm('N', 'N', this%mat_size, this%mat_size, this%mat_size, alfa, &
              bmat%matrix, bmat%mat_size, cmat%matrix, cmat%mat_size, beta,&
              this%matrix, this%mat_size)
#endif

   if (present(stat)) stat = my_stat
end subroutine mat_mul_x