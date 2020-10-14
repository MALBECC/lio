#include "datatypes/datatypes.fh"
module properties_data
   implicit none

   character(len=40) :: fmulliken
   character(len=40) :: fdipole
   character(len=40) :: flowdin
   character(len=40) :: fbecke

   logical :: fukui  = .false.
   logical :: dipole = .false.
   logical :: lowdin = .false.
   logical :: becke  = .false.

contains
end module properties_data