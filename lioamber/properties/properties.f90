#include "datatypes/datatypes.fh"
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
   public :: print_mulliken
   public :: print_becke
   public :: print_lowdin
   public :: print_fukui

contains

#include "printing_common.f90"
#include "mulliken.f90"
#include "lowdin.f90"
#include "becke.f90"
#include "misc.f90"
#include "fukui.f90"
end module properties