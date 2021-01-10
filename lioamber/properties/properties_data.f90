#include "../datatypes/datatypes.fh"
module properties_data
   implicit none

   character(len=40) :: fmulliken = "mulliken"
   character(len=40) :: fdipole   = "dipole_moment"
   character(len=40) :: flowdin   = "lowdin"
   character(len=40) :: fbecke    = "becke"
   character(len=40) :: ffukui    = "fukui"

   logical :: dipole   = .false.
   logical :: fukui    = .false.
   logical :: becke    = .false.
   logical :: lowdin   = .false.
   logical :: mulliken = .false.

   ! Internal UIDs for printing. For region printing,
   ! +50 is added to each UID.
   type uids_base
      integer :: dip  = 2000
      integer :: bec  = 2001
      integer :: low  = 2002
      integer :: mul  = 2003
      integer :: becs = 2004
      integer :: lows = 2005
      integer :: muls = 2006
      integer :: fuk  = 2007
      integer :: diptd= 2008
   end type uids_base
   type(uids_base) :: UIDs

   ! Data for properties grouped by region.
   type regions_base
      integer :: n_regions = 0           ! Number of regions.

      integer, allocatable :: natoms(:)  ! Number of atoms in a region.
      integer, allocatable :: atoms(:,:) ! Atom indexes for a given region.

      integer, allocatable :: nfuncs(:)  ! Number of basis functions in a region.
      integer, allocatable :: funcs(:,:) ! Basis indexes for a given region.
   end type regions_base
   type(regions_base) prop_regions
   
contains
end module properties_data