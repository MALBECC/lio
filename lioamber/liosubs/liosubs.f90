!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module liosubs
!------------------------------------------------------------------------------!
!
!  This is a place for general use subroutines that don't belong in more 
!  specific modules. It can also be used to temporarily store subroutines
!  that would later will packed in another module which is not ready yet.
!
!  Whatever the case, procedures here should not depend on other lio modules
!  but be as independent as possible instead.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! INTERFACES GO HERE (EXPLICIT)
   implicit none
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! PROCEDURES GO HERE (INCLUDED)
contains

#  include "catch_error.f90"
#  include "set_masses.f90"
#  include "nuclear_verlet.f90"

!  Possible IO-specific mode? 
#  include "catch_iostat.f90"
#  include "find_free_unit.f90"
#  include "safeio_open.f90"
#  include "safeio_rewind.f90"
#  include "write_geom.f90"
#  include "write_energy.f90"

end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
