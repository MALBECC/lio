!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module ehrendata
!------------------------------------------------------------------------------!
   implicit none
   integer :: nustep_count = 0
   integer :: elstep_count = 0

   real*8  :: StoredEnergy = 0.0d0
   integer :: rsti_funit  = 654321
   integer :: rsto_funit  = 123456

   complex*16,allocatable,dimension(:,:) :: RhoSaveA, RhoSaveB

end module ehrendata
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
