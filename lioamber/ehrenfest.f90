!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module ehrenfest
!------------------------------------------------------------------------------!
   implicit none
   real*8  :: StoredEnergy = 0.0d0
   integer :: step_number  = 0
   integer :: save_lapse   = 100
   integer :: last_step    = 1+120
   integer :: rstinp_unit  = 654321
   integer :: rstout_unit  = 123456
   logical :: save_step    = .false.
   logical :: restart_dyn  = .false.

   complex*16,allocatable,dimension(:,:) :: RhoSaveA, RhoSaveB
!
!
! INCLUDE FILES WITH HEADERS
!------------------------------------------------------------------------------!
#  include "ehrenfest/calc_Dmat_h.f90"
!
!
! INCLUDE FILES WITH PROCEDURES
!------------------------------------------------------------------------------!
contains

#  include "ehrenfest/ehrendyn.f90"
#  include "ehrenfest/ehrensetup.f90"

#  include "ehrenfest/setim.f90"
#  include "ehrenfest/calc_forceDS.f90"
#  include "ehrenfest/calc_forceDS_dss.f90"
#  include "ehrenfest/calc_forceDS_dds.f90"

#  include "ehrenfest/calc_Dmat_cholesky.f90"

#  include "ehrenfest/ehrenrst.f90"
#  include "ehrenfest/ehren_cholesky.f90"
#  include "ehrenfest/ehren_magnus.f90"
#  include "ehrenfest/ehren_verlet_e.f90"
#  include "ehrenfest/ehren_masses.f90"
#  include "ehrenfest/calc_kenergy.f90"

#  include "ehrenfest/RMMcalc0_Init.f90"
#  include "ehrenfest/RMMcalc1_Overlap.f90"
#  include "ehrenfest/RMMcalc2_FockMao.f90"

end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
