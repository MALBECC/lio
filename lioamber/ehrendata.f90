!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module ehrendata
!------------------------------------------------------------------------------!
   implicit none
   real*8  :: StoredEnergy = 0.0d0
   integer :: step_number  = 0
   integer :: rsti_funit  = 654321
   integer :: rsto_funit  = 123456

   complex*16,allocatable,dimension(:,:) :: RhoSaveA, RhoSaveB

end module ehrendata
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
