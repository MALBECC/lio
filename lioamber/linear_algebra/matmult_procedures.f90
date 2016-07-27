!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
function matmult_rr( mata, matb, option, info ) result( matc )
  implicit none
  real(kind=SPK), intent(in)     :: mata(:,:)
  real(kind=SPK), intent(in)     :: matb(:,:)
  real(kind=SPK)                 :: matc(size(mata,1),size(matb,2))
  integer, intent(in), optional  :: option
  integer, intent(out), optional :: info
# include "matmult_body.f90"
end function
!--------------------------------------------------------------------!
function matmult_dd( mata, matb, option, info ) result( matc )
  implicit none
  real(kind=DPK), intent(in)     :: mata(:,:)
  real(kind=DPK), intent(in)     :: matb(:,:)
  real(kind=DPK)                 :: matc(size(mata,1),size(matb,2))
  integer, intent(in), optional  :: option
  integer, intent(out), optional :: info
# include "matmult_body.f90"
end function
!--------------------------------------------------------------------!
function matmult_cc( mata, matb, option, info ) result( matc )
  implicit none
  complex(kind=SPK), intent(in)  :: mata(:,:)
  complex(kind=SPK), intent(in)  :: matb(:,:)
  complex(kind=SPK)              :: matc(size(mata,1),size(matb,2))
  integer, intent(in), optional  :: option
  integer, intent(out), optional :: info
# include "matmult_body.f90"
end function
!--------------------------------------------------------------------!
function matmult_zz( mata, matb, option, info ) result( matc )
  implicit none
  complex(kind=DPK), intent(in)  :: mata(:,:)
  complex(kind=DPK), intent(in)  :: matb(:,:)
  complex(kind=DPK)              :: matc(size(mata,1),size(matb,2))
  integer, intent(in), optional  :: option
  integer, intent(out), optional :: info
# include "matmult_body.f90"
end function
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
function matmult_rc( mata, matb, option, info ) result( matc )
  implicit none
  real(kind=SPK), intent(in)     :: mata(:,:)
  complex(kind=SPK), intent(in)  :: matb(:,:)
  complex(kind=SPK)              :: matc(size(mata,1),size(matb,2))
  integer, intent(in), optional  :: option
  integer, intent(out), optional :: info
# include "matmult_body.f90"
end function
!--------------------------------------------------------------------!
function matmult_cr( mata, matb, option, info ) result( matc )
  implicit none
  complex(kind=SPK), intent(in)  :: mata(:,:)
  real(kind=SPK), intent(in)     :: matb(:,:)
  complex(kind=SPK)              :: matc(size(mata,1),size(matb,2))
  integer, intent(in), optional  :: option
  integer, intent(out), optional :: info
# include "matmult_body.f90"
end function
!--------------------------------------------------------------------!
function matmult_dz( mata, matb, option, info ) result( matc )
  implicit none
  real(kind=DPK), intent(in)     :: mata(:,:)
  complex(kind=DPK), intent(in)  :: matb(:,:)
  complex(kind=DPK)              :: matc(size(mata,1),size(matb,2))
  integer, intent(in), optional  :: option
  integer, intent(out), optional :: info
# include "matmult_body.f90"
end function
!--------------------------------------------------------------------!
function matmult_zd( mata, matb, option, info ) result( matc )
  implicit none
  complex(kind=DPK), intent(in)  :: mata(:,:)
  real(kind=DPK), intent(in)     :: matb(:,:)
  complex(kind=DPK)              :: matc(size(mata,1),size(matb,2))
  integer, intent(in), optional  :: option
  integer, intent(out), optional :: info
# include "matmult_body.f90"
end function
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
function matmult_rd( mata, matb, option, info ) result( matc )
  real(kind=SPK), intent(in)     :: mata(:,:)
  real(kind=DPK), intent(in)     :: matb(:,:)
  real(kind=DPK)                 :: matc(size(mata,1),size(matb,2))
  integer, intent(in), optional  :: option
  integer, intent(out), optional :: info
# include "matmult_body.f90"
end function
!--------------------------------------------------------------------!
function matmult_dr( mata, matb, option, info ) result( matc )
  real(kind=DPK), intent(in)     :: mata(:,:)
  real(kind=SPK), intent(in)     :: matb(:,:)
  real(kind=DPK)                 :: matc(size(mata,1),size(matb,2))
  integer, intent(in), optional  :: option
  integer, intent(out), optional :: info
# include "matmult_body.f90"
end function
!--------------------------------------------------------------------!
function matmult_cz( mata, matb, option, info ) result( matc )
  implicit none
  complex(kind=SPK), intent(in)  :: mata(:,:)
  complex(kind=DPK), intent(in)  :: matb(:,:)
  complex(kind=DPK)              :: matc(size(mata,1),size(matb,2))
  integer, intent(in), optional  :: option
  integer, intent(out), optional :: info
# include "matmult_body.f90"
end function
!--------------------------------------------------------------------!
function matmult_zc( mata, matb, option, info ) result( matc )
  implicit none
  complex(kind=DPK), intent(in)  :: mata(:,:)
  complex(kind=SPK), intent(in)  :: matb(:,:)
  complex(kind=DPK)              :: matc(size(mata,1),size(matb,2))
  integer, intent(in), optional  :: option
  integer, intent(out), optional :: info
# include "matmult_body.f90"
end function
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
function matmult_rz( mata, matb, option, info ) result( matc )
  implicit none
  real(kind=SPK), intent(in)     :: mata(:,:)
  complex(kind=DPK), intent(in)  :: matb(:,:)
  complex(kind=DPK)              :: matc(size(mata,1),size(matb,2))
  integer, intent(in), optional  :: option
  integer, intent(out), optional :: info
# include "matmult_body.f90"
end function
!--------------------------------------------------------------------!
function matmult_zr( mata, matb, option, info ) result( matc )
  implicit none
  complex(kind=DPK), intent(in)  :: mata(:,:)
  real(kind=SPK), intent(in)     :: matb(:,:)
  complex(kind=DPK)              :: matc(size(mata,1),size(matb,2))
  integer, intent(in), optional  :: option
  integer, intent(out), optional :: info
# include "matmult_body.f90"
end function
!--------------------------------------------------------------------!
function matmult_dc( mata, matb, option, info ) result( matc )
  implicit none
  real(kind=DPK), intent(in)     :: mata(:,:)
  complex(kind=SPK), intent(in)  :: matb(:,:)
  complex(kind=DPK)              :: matc(size(mata,1),size(matb,2))
  integer, intent(in), optional  :: option
  integer, intent(out), optional :: info
# include "matmult_body.f90"
end function
!--------------------------------------------------------------------!
function matmult_cd( mata, matb, option, info ) result( matc )
  implicit none
  complex(kind=SPK), intent(in)  :: mata(:,:)
  real(kind=DPK), intent(in)     :: matb(:,:)
  complex(kind=DPK)              :: matc(size(mata,1),size(matb,2))
  integer, intent(in), optional  :: option
  integer, intent(out), optional :: info
# include "matmult_body.f90"
end function
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
