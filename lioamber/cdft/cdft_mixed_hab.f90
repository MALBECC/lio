
subroutine cdft_mixed_hab(e1, e2, Wat, Sat, is_open_shell)
   use cdft_data, only: cdft_mc, cdft_c, cdft_reg
   implicit none
   LIODBLE, intent(in)    :: e1, e2
   logical, intent(in)    :: is_open_shell
   LIODBLE, intent(inout) :: Wat(:,:), Sat(:,:)

   integer :: Msize, Nocup, Nocup2
   LIODBLE :: S_ab
   LIODBLE, allocatable :: Dmat(:,:)  , Dmat_inv(:,:), Omat(:,:)
   LIODBLE, allocatable :: Dmat_b(:,:), Dmat_inv_b(:,:), Omat_b(:,:)
   LIODBLE, allocatable :: tmpmat(:,:)

   Msize = size(Sat,1)
   Nocup = size(cdft_mc%coefs_a1,2)
   if (is_open_shell) Nocup2 = size(cdft_mc%coefs_b1,2)

   ! We first compute the easy term, Sab
   ! For this we first need the D matrix, which is Ca * S * Cb'
   call DGEMM('N', 'N', Msize, Nocup, Msize, 1.0D0, Sat, &
              Msize, cdft_mc%coefs_a2, Msize, 0.0D0, tmpmat, Msize)
   call DGEMM('N', 'T', Msize, Msize, Nocup, 1.0D0, cdft_mc%coefs_a1, &
              Msize, tmpmat, Msize, 0.0D0, Dmat, Msize)

   S_ab = get_determinant_qr(Dmat)

   if (is_open_shell) then
      call DGEMM('N', 'N', Msize, Nocup2, Msize, 1.0D0, Sat, &
                  Msize, cdft_mc%coefs_b2, Msize, 0.0D0, tmpmat, Msize)
      call DGEMM('N', 'T', Msize, Msize, Nocup2, 1.0D0, cdft_mc%coefs_b1, &
                  Msize, tmpmat, Msize, 0.0D0, Dmat_b, Msize)

      S_ab = S_ab * get_determinant_qr(Dmat_b)
   endif

   ! We now invert the matrix D.
   !call invert_matrix()

end subroutine cdft_mixed_hab

! Performs a QR decompostion to obtain the absolute value of the
! determinant of a matrix. Since the matrix we use is an overlap
! matrix, we do not really care about the sign. We also asume the
! matrix is square (spoiler alert: it is).
! In concrete:
!     S = Q*R (this is the decomposition)
! det(S) = det(Q)det(R), but Q is unitary so |det(Q)| = 1:
!     |det(S)| = |det(R)|
function get_determinant_qr(the_matrix) result(determinant)
   implicit none

   LIODBLE, intent(in) :: the_matrix(:,:)

   LIODBLE :: determinant
   integer :: mat_size, lwork, info, ii

   LIODBLE, allocatable :: bak_mat(:,:), tau(:), work(:)

   mat_size = size(the_matrix,1)

   ! We make a copy of the matrix, since DGEQRF destroys it.
   allocate(bak_mat(mat_size,mat_size))
   bak_mat = the_matrix

   ! Tau is an array needed to reconstruct Q.
   ! Since we do not care about Q, we leave it be.
   allocate(tau(mat_size))

   ! We guess optimal work size.
   allocate(work(1))
   call DGEQRF(mat_size, mat_size, bak_mat, mat_size, tau, &
               work, -1, info)

   lwork = int(work(1))
   deallocate(work); allocate(work(lwork))

   ! And finally, we perform the decomposition.
   call DGEQRF(mat_size, mat_size, bak_mat, mat_size, tau, &
               work, lwork, info)

   ! Now bak_mat contains R, which is upper triangular.
   ! This means the determinant is just the product of the
   ! diagonal elements.
   determinant = 1.0D0
   do ii = 1, mat_size
      determinant = determinant * bak_mat(ii,ii)
   enddo

   deallocate(bak_mat, tau, work)
end function get_determinant_qr