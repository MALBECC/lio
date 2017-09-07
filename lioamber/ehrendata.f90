!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module ehrendata
!------------------------------------------------------------------------------!
   implicit none
   integer :: nustep_count = 0
   integer :: elstep_count = 0

   integer :: rsti_funit  = 654321
   integer :: rsto_funit  = 123456

   real*8  :: stored_time   = 0.0d0
   real*8  :: stored_energy = 0.0d0
   real*8  :: stored_dipmom(3) = 0.0d0

   complex*16, allocatable :: stored_densM1(:,:)
   complex*16, allocatable :: stored_densM2(:,:)

end module ehrendata
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
