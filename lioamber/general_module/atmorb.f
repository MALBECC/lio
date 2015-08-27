!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine atmorb_n(atmvec,atm_of_orb,orbvec)
       implicit none
       integer,intent(in)   :: atmvec(:)
       integer,intent(in)   :: atm_of_orb(:)
       integer,intent(out)  :: orbvec(:)
       integer :: kk

!------------------------------------------------------------------------------!
! CAN'T DO CHECK BECAUSE OF LEGACY PRACTICE OF SETTING FIXED SIZES
!       if (size(atm_of_orb).ne.size(orbvec))
!     >   stop ('atmorb: Problem with size of basis set.')
!       if (maxval(atm_of_orb).ne.size(atmvec))
!     >   stop ('atmorb: Problem with size of atoms set.')

       do kk=1,size(orbvec)
         orbvec(kk)=atmvec(atm_of_orb(kk))
       enddo

!------------------------------------------------------------------------------!
       return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
