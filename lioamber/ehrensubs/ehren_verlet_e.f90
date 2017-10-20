!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine ehren_verlet_e(Nsize,dt,Fmat,Rold,Rnow,Rnew)
!
! Fmat,Rold,Rnew => All in the ON basis
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in)     :: Nsize
  real*8,intent(in)      :: dt
!  real*8,intent(in)      :: Fmat(Nsize,Nsize)
  complex*16,intent(in)  :: Fmat(Nsize,Nsize)
  complex*16,intent(in)  :: Rold(Nsize,Nsize)
  complex*16,intent(in)  :: Rnow(Nsize,Nsize)
  complex*16,intent(out) :: Rnew(Nsize,Nsize)

  complex*16,allocatable :: ConmMat(:,:)
  complex*16,allocatable :: TermPos(:,:)
  complex*16,allocatable :: TermNeg(:,:)


  allocate(ConmMat(Nsize,Nsize),TermPos(Nsize,Nsize),TermNeg(Nsize,Nsize))
  TermPos=matmul(Fmat,Rnow)
  TermNeg=matmul(Rnow,Fmat)
  ConmMat=TermPos-TermNeg
  Rnew=Rold-dt*CMPLX(0.0d0,2.0d0)*ConmMat
!  ConmMat=ConmMat*CMPLX(0.0d0,-1.0d0)
!  Rnew=Rold+(2.0d0)*dt*ConmMat

  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
