!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module converger_data

   implicit none

!  Covergens criterion
!  Damping=1, DIIS=2, Hybrid Converg=3
   integer :: conver_criter

   logical :: hagodiis
   integer :: ndiis
   real*8  :: damping_factor

   real*8, allocatable :: fock_damped(:,:,:)
   real*8, allocatable :: bcoef (:,:)
   real*8, allocatable :: fockm (:,:,:,:)
   real*8, allocatable :: FP_PFm (:,:,:,:)
   real*8, allocatable :: EMAT2 (:,:,:)

end module converger_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
