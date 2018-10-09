!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module liosubs
!------------------------------------------------------------------------------!
!
!  This is a place for general use subroutines that don't belong in more 
!  specific modules. It can also be used to temporarily store subroutines
!  that would later will packed in another module which is not ready yet.
!  Whatever the case, procedures here should not depend on other lio modules
!  but be as independent as possible instead.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! INTERFACES GO HERE (EXPLICIT)
   implicit none

   interface check_vecsize
      module procedure check_vecsize_n
      module procedure check_vecsize_r
      module procedure check_vecsize_d
      module procedure check_vecsize_c
      module procedure check_vecsize_z
   end interface check_vecsize

   interface check_matsize
      module procedure check_matsize_n
      module procedure check_matsize_r
      module procedure check_matsize_d
      module procedure check_matsize_c
      module procedure check_matsize_z
   end interface check_matsize

   interface atmvec_to_orbvec
      module procedure atmvec_to_orbvec_n
      module procedure atmvec_to_orbvec_r
      module procedure atmvec_to_orbvec_d
      module procedure atmvec_to_orbvec_c
      module procedure atmvec_to_orbvec_z
   end interface atmvec_to_orbvec

   interface read_list
      module procedure read_list_n
      module procedure read_list_r
      module procedure read_list_d
      module procedure read_list_c
      module procedure read_list_z
   end interface read_list

!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! PROCEDURES GO HERE (INCLUDED)
contains

#  include "line_search.f90"
#  include "catch_error.f90"
#  include "check_vecsize.header.f90"
#  include "check_matsize.header.f90"
#  include "atmvec_to_orbvec.header.f90"
#  include "read_list.header.f90"

#  include "gaussian_shaper.f90"
#  include "set_masses.f90"
#  include "nuclear_verlet.f90"

!  Possible IO-specific module? Think if the following procedures shouldn't
!  go together but somewhere else.
#  include "catch_iostat.f90"
#  include "find_free_unit.f90"
#  include "safeio_open.f90"
#  include "safeio_rewind.f90"
#  include "write_geom.f90"
#  include "write_energy.f90"

end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
