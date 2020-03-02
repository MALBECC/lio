#include "datatypes/datatypes.fh"
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module fockbias_data
!------------------------------------------------------------------------------!
!
! DESCRIPTION PENDING
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none

!  Variables that are read from input or forced init (only once)
   logical             :: fockbias_is_active = .false.
   logical             :: fockbias_is_shaped = .false.
   LIODBLE              :: fockbias_timegrow  = 0.0d0
   LIODBLE              :: fockbias_timefall  = 0.0d0
   LIODBLE              :: fockbias_timeamp0  = 0.0d0
   character(len=80)   :: fockbias_readfile  = "atombias.in"

!  Variables that are configured during setorb (only once)
   LIODBLE, allocatable :: fockbias_orbqw(:)

!  Variables that are configured during setmat (redo if atoms move)
   LIODBLE, allocatable :: fockbias_matrix(:,:)

end module fockbias_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
