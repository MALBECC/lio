!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!subroutine atmvec_to_orbvec_X( atmvec, atm_of_orb, orbvec )
!   implicit none
!   XXXXXXX,intent(in)  :: atmvec(:)
!   integer,intent(in)  :: atm_of_orb(:)
!   XXXXXXX,intent(out) :: orbvec(:)
!------------------------------------------------------------------------------!
   integer :: nn

   nn = size(atm_of_orb)
   call check_vecsize( nn, orbvec, "orbvec", "atmvec_to_orbvec" )

   nn = max( size(atmvec), maxval(atm_of_orb) )
   call check_vecsize( nn, atmvec, "atmvec", "atmvec_to_orbvec" )

   do nn = 1, size(orbvec)
      orbvec(nn) = atmvec( atm_of_orb(nn) )
   enddo

!------------------------------------------------------------------------------!
!end subroutine atmvec_to_orbvec_X
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
