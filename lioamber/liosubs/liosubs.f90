!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module liosubs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
! This is a place for general use subroutines that don't belong
! in more specific modules. Procedures here should not depend on
! other lio modules but be as independent as possible instead.
!
!--------------------------------------------------------------------!
  implicit none
  contains

#   include "set_masses.f90"
#   include "nuclear_verlet.f90"

! The following procedures should be exported to a io-specific mod
#   include "catch_iostat.f90"
#   include "find_free_unit.f90"
#   include "write_geom.f90"
#   include "write_energy.f90"

end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
