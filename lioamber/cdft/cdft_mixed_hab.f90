! Calculates orthogonalised Hab values.
! We first obtain the non-orthogonalised Hab values from the following
! equivalence:
!
! Hab = 0.5 * (Fa + Fb) * Sab + SUM_{regions} 0.5 * < a | Wr_a + Wr_b | b >
!
! We define 0.5 * < a | Wr_a + Wr_b | b > as the W matrix.
!
! This is, in truth, the average of Hab and Hba, but it guarantees that
! the H matrix in the {a,b} is Hermitian. We then orthogonalise the H
! matrix by transforming H via V* H V, where V is the eigenvector matrix
! of the W matrix. In concrete, we obtain V by solving:
!     W * V = S * V * E
! Where E is the eigenvalue diagonal matrix.
subroutine cdft_mixed_hab(Ea, Eb, Wat, Wat_b, Sat, is_open_shell, Hmat, S_ab)
   use cdft_data, only: cdft_mc

   implicit none
   LIODBLE, intent(in)    :: Ea, Eb
   logical, intent(in)    :: is_open_shell
   LIODBLE, intent(inout) :: Wat(:,:), Wat_b(:,:), Sat(:,:), Hmat(2,2), S_ab

   integer :: Msize, Nocup, Nocup2
   integer :: iorb, jorb
   LIODBLE :: H_ab, accum_o
   LIODBLE, allocatable :: Dmat(:,:)  , Dmat_inv(:,:), Omat(:,:)
   LIODBLE, allocatable :: Dmat_b(:,:), Dmat_inv_b(:,:), Omat_b(:,:)
   LIODBLE, allocatable :: tmpmat(:,:)

   Msize = size(Sat,1)
   Nocup = size(cdft_mc%coefs_a1,2)
   if (is_open_shell) Nocup2 = size(cdft_mc%coefs_b1,2)

   ! We first compute the easy term, Sab
   ! For this we first need the D matrix, which is Ca * S * Cb'
   allocate(Dmat(Nocup,Nocup))
   allocate(tmpmat(Msize,Nocup))
   tmpmat = 0.0D0
   Dmat   = 0.0D0
   call DGEMM('N', 'N', Msize, Nocup, Msize, 1.0D0, Sat, &
              Msize, cdft_mc%coefs_a1, Msize, 0.0D0, tmpmat, Msize)
   call DGEMM('T', 'N', Nocup, Nocup, Msize, 1.0D0, cdft_mc%coefs_a2, &
              Msize, tmpmat, Msize, 0.0D0, Dmat, Nocup)
   deallocate(tmpmat)

   S_ab = get_determinant_qr(Dmat)

   ! We alloc-dealloc-realloc open shell variables to avoid
   ! compilation warnings.
   allocate(Dmat_b(1,1))
   if (is_open_shell) then
      allocate(tmpmat(Msize,Nocup2))
      deallocate(Dmat_b); allocate(Dmat_b(Nocup2,Nocup2))
      tmpmat = 0.0D0
      Dmat_b = 0.0D0
   
      call DGEMM('N', 'N', Msize, Nocup2, Msize, 1.0D0, Sat, &
                  Msize, cdft_mc%coefs_b1, Msize, 0.0D0, tmpmat, Msize)
      call DGEMM('T', 'N', Nocup2, Nocup2, Msize, 1.0D0, cdft_mc%coefs_b2, &
                  Msize, tmpmat, Msize, 0.0D0, Dmat_b, Nocup2)
      deallocate(tmpmat)

      accum_o = get_determinant_qr(Dmat_b)
      S_ab = S_ab * accum_o
   endif

   ! We now invert the matrix D. Afterwards, we do not need D any longer.
   allocate(Dmat_inv(Nocup,Nocup))
   Dmat_inv = 0.0D0
   call get_inverse_matrix(Dmat, Dmat_inv)
   deallocate(Dmat)

   allocate(Dmat_inv_b(1,1))
   if (is_open_shell) then
      deallocate(Dmat_inv_b); allocate(Dmat_inv_b(Nocup2,Nocup2))
      Dmat_inv_b = 0.0D0
      call get_inverse_matrix(Dmat_b, Dmat_inv_b)
      deallocate(Dmat_b)
   endif

   ! We now obtain our last ingredient: the O matrix, which is Ca * W * Cb'
   allocate(Omat(Nocup,Nocup))
   allocate(tmpmat(Msize,Nocup))
   Omat   = 0.0D0
   tmpmat = 0.0D0

   call DGEMM('N', 'N', Msize, Nocup, Msize, 1.0D0, Wat, &
               Msize, cdft_mc%coefs_a1, Msize, 0.0D0, tmpmat, Msize)
   call DGEMM('T', 'N', Nocup, Nocup, Msize, 1.0D0, cdft_mc%coefs_a2, &
               Msize, tmpmat, Msize, 0.0D0, Omat, Nocup)
   deallocate(tmpmat)

   allocate(Omat_b(1,1))
   if (is_open_shell) then
      deallocate(Omat_b); allocate(Omat_b(Nocup2,Nocup2))
      allocate(tmpmat(Msize,Nocup2))
      Omat_b = 0.0D0
      tmpmat = 0.0D0
      call DGEMM('N', 'N', Msize, Nocup2, Msize, 1.0D0, Wat_b, &
                 Msize, cdft_mc%coefs_b1, Msize, 0.0D0, tmpmat, Msize)
      call DGEMM('T', 'N', Nocup2, Nocup2, Msize, 1.0D0, cdft_mc%coefs_b1, &
                 Msize, tmpmat, Msize, 0.0D0, Omat_b, Nocup2)
      deallocate(tmpmat)
   endif

   ! And now we can finally achieve greatness.
   H_ab = 0.5D0 * (Ea + Eb) * S_ab

   accum_o = 0.0D0
   do iorb = 1, Nocup
   do jorb = 1, Nocup
      accum_o = accum_o + abs(Omat(iorb,jorb) * Dmat_inv(iorb,jorb))
   enddo
   enddo
   
   if (is_open_shell) then
      do iorb = 1, Nocup2
      do jorb = 1, Nocup2
         accum_o = accum_o + Omat_b(iorb,jorb) * Dmat_inv_b(iorb,jorb)
      enddo
      enddo
   endif
   H_ab = H_ab - S_ab * accum_o

   ! But this is not even my final form! We still need to obtain the proper
   ! projection, so that we find the Hamiltonian in an orthogonal basis.
   ! We must perform the following transformation (Lowdin orthogonalisation):
   ! H' = S^-1/2 * H * S^-1/2
   call orthog_Hab(H_ab, S_ab, Ea, Eb, accum_o, Hmat)
   deallocate(Dmat_inv, Omat)
   deallocate(Dmat_inv_b, Omat_b)
end subroutine cdft_mixed_hab

! Orthogonolises H matrix in the basis {A,B} by doing
! H' = S^-1/2 * H * S^-1/2
subroutine orthog_Hab(Hab, Sab, Ea, Eb, Wab, Hortho)
   use cdft_data, only: cdft_reg, cdft_c
   implicit none
   LIODBLE, intent(in)  :: Sab, Ea, Eb, Hab, Wab
   LIODBLE, intent(out) :: Hortho(2,2)

   logical :: easy_orthog  = .true.
   LIODBLE, allocatable :: Smat(:,:), Hmat(:,:)
   LIODBLE, allocatable :: Umat(:,:), tmpmat(:,:)
   LIODBLE, allocatable :: evals(:), work(:)
   integer :: info, ii

   ! Yet another orthogonalisation. We do W*C = S*C*n, where
   ! n is a diagonal matrix containing the eigenvalues. 
   ! We use DSYGV to solve the generalised eigenvalue problem.

   allocate(Smat(2,2)) ! S matrix
   Smat = reshape((/1.0D0, Sab, Sab, 1.0D0/), shape(Smat))
   
   allocate(tmpmat(2,2), umat(2,2))
   umat = 0.0D0
   do ii = 1, cdft_c%n_regions
      umat(1,1) = umat(1,1) + cdft_reg%Vc2(ii) * &
                              (cdft_reg%nelecs(ii) - cdft_reg%chrg2(ii))
      umat(1,1) = umat(1,1) + cdft_reg%Vc(ii) * &
                              (cdft_reg%nelecs(ii) - cdft_reg%chrg(ii))
   enddo
   umat(2,1) = Wab; umat(1,2) = Wab

   tmpmat = umat

   allocate(evals(2), work(10))
   call dsygv(1, 'V', 'U', 2, tmpmat, 2, Smat, 2, evals, work, 10, info)
   deallocate(evals, work)

   allocate(Hmat(2,2))

   ! We obtain H' as shown in 10.10163/1.3507878
   Hmat = reshape((/Ea, Hab, Hab, Eb/), shape(Hmat))
   
   Hortho = matmul(Hmat,tmpmat)
   tmpmat = transpose(tmpmat)
   Hortho = matmul(tmpmat, Hortho)

   deallocate(Smat, tmpmat, umat, Hmat)
   return

   ! Silly orthogonalisation as shown in DOI:10.1039/C7CP06660K
   ! I will keep this here for debugging purposes.
   Hortho = 0.0D0
   if (easy_orthog) then
      Hortho(1,1) = Ea
      Hortho(2,2) = Eb

      Hortho(1,2) = (Hab - 0.5D0 * Sab * (Ea + Eb)) / (1.0D0 - Sab * Sab)
      Hortho(2,1) = Hortho(1,2)
      return
   endif
 end subroutine orthog_Hab

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
   if (info /= 0) print*, "Problem in determinant. Error = ", info

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
   determinant = abs(determinant)

   deallocate(bak_mat, tau, work)
end function get_determinant_qr

! Performs the inverse of a matrix and stores the output in
! out_matrix.
subroutine get_inverse_matrix(the_matrix, out_matrix)
   implicit none
   ! DGETRF y luego DGETRI
   LIODBLE, intent(in)  :: the_matrix(:,:)
   LIODBLE, intent(out) :: out_matrix(:,:)

   integer :: mat_size, info, lwork
   
   integer, allocatable :: pivot_index(:)
   LIODBLE, allocatable :: work(:)

   mat_size = size(the_matrix,1)

   allocate(pivot_index(mat_size))
   ! We make a copy of the matrix, since DGETRF destroys it.
   out_matrix = the_matrix

   ! Make P*L*U factorization.
   call DGETRF(mat_size, mat_size, out_matrix, mat_size, pivot_index, info)

   ! Calculate inverse using the PLU factorization.
   allocate(work(1))
   call DGETRI(mat_size, out_matrix, mat_size, pivot_index, work, -1, info)

   lwork = int(work(1))
   deallocate(work); allocate(work(lwork))
   call DGETRI(mat_size, out_matrix, mat_size, pivot_index, work, lwork, info)
   deallocate(work)
end subroutine get_inverse_matrix
