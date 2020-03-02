!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrenaux_cholesky( Nsize, Smat, Lmat, Umat, Linv, Uinv, Sinv )
!------------------------------------------------------------------------------!
!
! DESCRIPTION
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in) :: Nsize
  LIODBLE,intent(in)  :: Smat(Nsize,Nsize)
  LIODBLE,intent(out) :: Lmat(Nsize,Nsize),Umat(Nsize,Nsize)
  LIODBLE,intent(out) :: Linv(Nsize,Nsize),Uinv(Nsize,Nsize)
  LIODBLE,intent(out) :: Sinv(Nsize,Nsize)

  integer :: ii,jj,iost


  Lmat=Smat
  call dpotrf('L', Nsize, Lmat, Nsize, iost)
  do ii=1,Nsize-1
  do jj=ii+1,Nsize
    Lmat(ii,jj)=0.0d0
  enddo
  enddo


  Linv=Lmat
  call dtrtri('L', 'N', Nsize, Linv, Nsize, iost)
  Umat=transpose(Lmat)
  Uinv=transpose(Linv)
  Sinv=matmul(Uinv,Linv)


end subroutine ehrenaux_cholesky
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
