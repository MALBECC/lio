module kinetic_data
  implicit none

  logical                   :: kin_separation =.true. ! Separtes the kinetic energy of the contribution of 1 e integrals
  real(kind=8), allocatable :: Tmat(:)                ! Kinetic energy matrix
  real(kind=8)              :: KinE                   ! Kinetic energy accumulator

end module kinetic_data
