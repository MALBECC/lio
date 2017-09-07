!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module ehrensubs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
#  include "calc_Dmat_h.f90"
!
!
contains
#  include "ehrendyn_main.f90"
#  include "ehrendyn_init.f90"
#  include "ehrendyn_step.f90"
#  include "ehrendyn_prep.f90"

#  include "calc_gintmat.f90"
#  include "calc_forceDS.f90"
#  include "calc_forceDS_dss.f90"
#  include "calc_forceDS_dds.f90"
#  include "calc_Dmat_cholesky.f90"

#  include "ehrenaux_rstx.f90"
#  include "ehrenaux_setfld.f90"
#  include "ehrenaux_cholesky.f90"
#  include "ehrenaux_magnus.f90"
#  include "ehrenaux_verlet.f90"
#  include "ehrenaux_masses.f90"
#  include "ehrenaux_calckyn.f90"
#  include "ehrenaux_writedip.f90"
#  include "ehrenaux_updatevel.f90"

#  include "RMMcalc0_Init.f90"
#  include "RMMcalc1_Overlap.f90"
#  include "RMMcalc2_FockMao.f90"
#  include "RMMcalc3_FockMao.f90"
#  include "RMMcalc4_FockMao.f90"

end module ehrensubs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
