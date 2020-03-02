!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrenaux_verlet( Nsize, dt, Fmat, Rold, Rnow, Rnew )
!------------------------------------------------------------------------------!
!
! DESCRIPTION
! Fmat,Rold,Rnew => All in the ON basis
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in)     :: Nsize
  LIODBLE,intent(in)      :: dt
  complex(kind=8),intent(in)  :: Fmat(Nsize,Nsize)
  complex(kind=8),intent(in)  :: Rold(Nsize,Nsize)
  complex(kind=8),intent(in)  :: Rnow(Nsize,Nsize)
  complex(kind=8),intent(out) :: Rnew(Nsize,Nsize)

  complex(kind=8),allocatable :: ConmMat(:,:)
  complex(kind=8),allocatable :: TermPos(:,:)
  complex(kind=8),allocatable :: TermNeg(:,:)


  allocate( ConmMat(Nsize,Nsize), TermPos(Nsize,Nsize), TermNeg(Nsize,Nsize) )
  TermPos = matmul(Fmat, Rnow)
  TermNeg = matmul(Rnow, Fmat)
  ConmMat = TermPos - TermNeg
  Rnew = Rold - dt * CMPLX(0.0d0,2.0d0,8) * ConmMat

end subroutine ehrenaux_verlet
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
