#include "datatypes/datatypes.fh"
module properties
   implicit none


   interface mulliken
      module procedure mulliken_cs
      module procedure mulliken_os
   end interface mulliken
contains

#include "mulliken.f90"
#include "lowdin.f90"
#include "misc.f90"
#include "fukui.f90"
end module properties