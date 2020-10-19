#include "../datatypes/datatypes.fh"
module properties
   implicit none
   private

   interface mulliken
      module procedure mulliken_cs
      module procedure mulliken_os
   end interface mulliken

   interface print_mulliken
      module procedure print_mulliken_cs
      module procedure print_mulliken_os
   end interface print_mulliken

   interface lowdin
      module procedure lowdin_cs
      module procedure lowdin_os
   end interface lowdin

   interface print_lowdin
      module procedure print_lowdin_cs
      module procedure print_lowdin_os
   end interface print_lowdin

   interface fukui
      module procedure fukui_calc_cs
      module procedure fukui_calc_os
   end interface fukui

   interface print_fukui
      module procedure print_fukui_cs
      module procedure print_fukui_os
   end interface print_fukui

   public :: lowdin
   public :: mulliken
   public :: fukui
   public :: do_lowdin
   public :: do_mulliken
   public :: do_becke
   public :: do_fukui
   public :: do_dipole
   public :: print_mulliken
   public :: print_becke
   public :: print_lowdin
   public :: print_fukui
   public :: dipole
   public :: write_dipole_td_header
   public :: write_dipole_td
   public :: properties_initialise
   public :: properties_finalise
   public :: properties_region_read
contains

#include "write_population.f90"
#include "mulliken.f90"
#include "lowdin.f90"
#include "becke.f90"
#include "misc.f90"
#include "fukui.f90"
#include "dipole.f90"
#include "init_fin.f90"
#include "regions_read.f90"


function do_mulliken() result(do_out)
   use properties_data, only: mulliken
   implicit none
   logical :: do_out
   do_out = mulliken
end function do_mulliken

function do_lowdin() result(do_out)
   use properties_data, only: lowdin
   implicit none
   logical :: do_out
   do_out = lowdin
end function do_lowdin

function do_becke() result(do_out)
   use properties_data, only: becke
   implicit none
   logical :: do_out
   do_out = becke
end function do_becke

function do_fukui() result(do_out)
   use properties_data, only: fukui
   implicit none
   logical :: do_out
   do_out = fukui
end function do_fukui

function do_dipole() result(do_out)
   use properties_data, only: dipole
   implicit none
   logical :: do_out
   do_out = dipole
end function do_dipole


end module properties