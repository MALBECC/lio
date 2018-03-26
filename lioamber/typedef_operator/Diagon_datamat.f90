subroutine Diagon_datamat (this, eigen_vecs, eigen_vals)
   use linear_algebra, only: matrix_diagon

   implicit none
   class(operator), intent(in) :: this
   real*8, intent(out) :: eigen_vecs(:,:)
   real*8, intent(out) :: eigen_vals(:)
   real*8, allocatable :: Dmat(:,:)
   integer             :: Nbasis

   Nbasis=size(this%data_ON,1)

   allocate(Dmat(Nbasis,Nbasis))

   Dmat=  this%data_ON 

   call matrix_diagon( Dmat, eigen_vecs, eigen_vals )

end subroutine Diagon_datamat
