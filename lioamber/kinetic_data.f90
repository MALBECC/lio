#include "datatypes/datatypes.fh"

module kinetic_data
  implicit none

  logical              :: kin_separation =.true.      ! Separtes the kinetic energy of the contribution of 1 e integrals
  LIODBLE, allocatable :: Tmat(:)                     ! Kinetic energy matrix
  LIODBLE, allocatable :: pro(:)                      ! Basis function i projection over density
  LIODBLE, allocatable :: Smatd(:,:)                  ! Auxiliary basis set overlap matrix
  LIODBLE              :: KinE                        ! Kinetic energy accumulator

end module kinetic_data
