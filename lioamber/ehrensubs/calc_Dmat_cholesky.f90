!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
function calc_Dmat_cholesky( nbasis, Lmat, Umat, Bmat ) result(Dmat)
!------------------------------------------------------------------------------!
!
! DESCRIPTION
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in) :: nbasis
  LIODBLE,intent(in)  :: Lmat(nbasis,nbasis) ! From Cholesky Decomp
  LIODBLE,intent(in)  :: Umat(nbasis,nbasis) ! From Cholesky Decomp
  LIODBLE,intent(in)  :: Bmat(nbasis,nbasis)
  LIODBLE             :: Dmat(nbasis,nbasis)

  LIODBLE,allocatable :: Matrix(:,:)
  integer            :: ii,jj

!------------------------------------------------------------------------------!

  allocate(Matrix(nbasis,nbasis))
  Matrix=(-1.0d0)*Lmat
  Matrix=matmul(Matrix,Bmat)
  Matrix=matmul(Matrix,Umat)


  do jj=1,nbasis

     do ii=1,jj-1
       Dmat(ii,jj)=-Matrix(jj,ii)
     enddo

     Dmat(jj,jj)=0.0d0

     do ii=jj+1,nbasis
       Dmat(ii,jj)=Matrix(ii,jj)
     enddo

  enddo


  deallocate(Matrix)
end function calc_Dmat_cholesky
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
