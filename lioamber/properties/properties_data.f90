#include "datatypes/datatypes.fh"
module properties_data
   implicit none

   character(len=40) :: fmulliken = "mulliken"
   character(len=40) :: fdipole   = "dipole"
   character(len=40) :: flowdin   = "lowdin"
   character(len=40) :: fbecke    = "becke"

   logical :: dipole   = .false.
   logical :: fukui    = .false.
   logical :: becke    = .false.
   logical :: lowdin   = .false.
   logical :: mulliken = .false.

   ! Internal UIDs for printing.
   type uids_base
      integer :: dip  = 2000
      integer :: bec  = 2001
      integer :: low  = 2002
      integer :: mul  = 2003
      integer :: becs = 2004
      integer :: lows = 2005
      integer :: muls = 2006
   end type uids_base
   type(uids_base) :: UIDs
   
contains
end module properties_data