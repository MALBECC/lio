!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
#include "../datatypes/datatypes.fh"
module ehrensubs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
#  include "calc_Dmat_h.f90"

interface  rmmget_dens
   module procedure rmmget_dens_r
   module procedure rmmget_dens_d
   module procedure rmmget_dens_c
   module procedure rmmget_dens_z
end interface rmmget_dens

interface  rmmput_dens
   module procedure rmmput_dens_r
   module procedure rmmput_dens_d
   module procedure rmmput_dens_c
   module procedure rmmput_dens_z
end interface rmmput_dens

interface  rmmget_fock
   module procedure rmmget_fock_r
   module procedure rmmget_fock_d
end interface rmmget_fock

interface  rmmput_fock
   module procedure rmmput_fock_r
   module procedure rmmput_fock_d
end interface rmmput_fock

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

#  include "rmmget_dens.f90"
#  include "rmmget_fock.f90"
#  include "rmmput_dens.f90"
#  include "rmmput_fock.f90"

#  include "RMMcalc0_Init.f90"
#  include "RMMcalc1_Overlap.f90"
#  include "RMMcalc2_FockMao.f90"
#  include "RMMcalc3_FockMao.f90"
#  include "RMMcalc4_FockMao.f90"

end module ehrensubs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
