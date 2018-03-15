subroutine Diagon_datamat (this, eigen_vecs, eigen_vals)
   use linear_algebra, only: matrix_diagon

   implicit none
   class(operator), intent(in) :: this
   real*8, intent(out) :: eigen_vecs(:,:)
   real*8, intent(out) :: eigen_vals(:)

   call matrix_diagon( this%data_ON, eigen_vecs, eigen_vals )

end subroutine Diagon_datamat
