!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrenaux_magnus( Nsize, Norder, dt, Fmat, Rold, Rnew )
!------------------------------------------------------------------------------!
!
! DESCRIPTION
! Fmat,Rold,Rnew => All in the ON basis
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in)     :: Nsize,Norder
  LIODBLE,intent(in)      :: dt
  complex(kind=8),intent(in)  :: Fmat(Nsize,Nsize)
  complex(kind=8),intent(in)  :: Rold(Nsize,Nsize)
  complex(kind=8),intent(out) :: Rnew(Nsize,Nsize)

  integer :: kk
  LIODBLE  :: factinv
  complex(kind=8),allocatable :: Omega1(:,:)
  complex(kind=8),allocatable :: ConmutAcum(:,:)
  complex(kind=8),allocatable :: TermPos(:,:),TermNeg(:,:)
!
!
! Initializations
!------------------------------------------------------------------------------!
  allocate( Omega1(Nsize,Nsize) )
  allocate( ConmutAcum(Nsize,Nsize) )
  allocate( TermPos(Nsize,Nsize), TermNeg(Nsize,Nsize) )

  Rnew = Rold
  ConmutAcum = Rold
  Omega1 = DCMPLX(0.0d0,-1.0d0) * (Fmat) * (dt)
!  factinv = 1.0d0
!
!
! Calculations
!------------------------------------------------------------------------------!
  do kk = 1, Norder
    TermPos = matmul(Omega1, ConmutAcum)
    TermNeg = matmul(ConmutAcum, Omega1)
    factinv = 1.0d0 / dble(kk)
    ConmutAcum = (TermPos - TermNeg) * factinv
    Rnew = Rnew + ConmutAcum
!    factinv = factinv / dble(kk)
!    Rnew = Rnew + factinv * ConmutAcum
  enddo

  deallocate( Omega1, ConmutAcum, TermPos, TermNeg )
end subroutine ehrenaux_magnus
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
