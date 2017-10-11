!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine atmvec_to_orbvec_n( atmvec, atm_of_orb, orbvec )
   implicit none
   integer   , intent(in)  :: atmvec(:)
   integer   , intent(in)  :: atm_of_orb(:)
   integer   , intent(out) :: orbvec(:)
#  include "atmvec_to_orbvec.proced.f90"
end subroutine atmvec_to_orbvec_n
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine atmvec_to_orbvec_r( atmvec, atm_of_orb, orbvec )
   implicit none
   real*4    , intent(in)  :: atmvec(:)
   integer   , intent(in)  :: atm_of_orb(:)
   real*4    , intent(out) :: orbvec(:)
#  include "atmvec_to_orbvec.proced.f90"
end subroutine atmvec_to_orbvec_r
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine atmvec_to_orbvec_d( atmvec, atm_of_orb, orbvec )
   implicit none
   real*8    , intent(in)  :: atmvec(:)
   integer   , intent(in)  :: atm_of_orb(:)
   real*8    , intent(out) :: orbvec(:)
#  include "atmvec_to_orbvec.proced.f90"
end subroutine atmvec_to_orbvec_d
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine atmvec_to_orbvec_c( atmvec, atm_of_orb, orbvec )
   implicit none
   complex*8 , intent(in)  :: atmvec(:)
   integer   , intent(in)  :: atm_of_orb(:)
   complex*8 , intent(out) :: orbvec(:)
#  include "atmvec_to_orbvec.proced.f90"
end subroutine atmvec_to_orbvec_c
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine atmvec_to_orbvec_z( atmvec, atm_of_orb, orbvec )
   implicit none
   complex*16, intent(in)  :: atmvec(:)
   integer   , intent(in)  :: atm_of_orb(:)
   complex*16, intent(out) :: orbvec(:)
#  include "atmvec_to_orbvec.proced.f90"
end subroutine atmvec_to_orbvec_z
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
