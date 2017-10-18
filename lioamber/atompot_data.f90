!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module atompot_data
!------------------------------------------------------------------------------!
!
! DESCRIPTION PENDING
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  logical             :: atompot_apply = .false.
  logical             :: atompot_ready = .false.
  integer             :: atompot_msize = 0
  real*8, allocatable :: qweight_of_orb(:)

  logical             :: atompot_timegrow = .false.
  logical             :: atompot_timefall = .false.
  real*8              :: atompot_timepos0 = 0.0d0
  real*8              :: atompot_timeamp1 = 0.0d0

end module atompot_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
