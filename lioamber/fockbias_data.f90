!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module fockbias_data
!------------------------------------------------------------------------------!
!
! DESCRIPTION PENDING
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  logical             :: fockbias_apply = .false.
  logical             :: fockbias_ready = .false.
  integer             :: fockbias_msize = 0
  real*8, allocatable :: qweight_of_orb(:)

  logical             :: fockbias_timegrow = .false.
  logical             :: fockbias_timefall = .false.
  real*8              :: fockbias_timepos0 = 0.0d0
  real*8              :: fockbias_timeamp1 = 0.0d0

end module fockbias_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
