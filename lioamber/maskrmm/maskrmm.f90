!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module maskrmm
   implicit none

!------------------------------------------------------------------------------!
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


!------------------------------------------------------------------------------!
   interface  rmmget_fock
      module procedure rmmget_fock_r
      module procedure rmmget_fock_d
   end interface rmmget_fock

   interface  rmmput_fock
      module procedure rmmput_fock_r
      module procedure rmmput_fock_d
   end interface rmmput_fock


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   contains
#     include "rmmput_dens.f90"
#     include "rmmget_dens.f90"
#     include "rmmput_fock.f90"
#     include "rmmget_fock.f90"

#     include "rmmcalc0_init.f90"
#     include "rmmcalc1_overlap.f90"
#     include "rmmcalc2_focknuc.f90"
#     include "rmmcalc3_fockele.f90"

end module maskrmm
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
