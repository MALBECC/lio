!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine catch_error( caller_name, statid_chck, statid_pass, statid_comm )
!------------------------------------------------------------------------------!
!
!     The following subroutine provides a simple, direct and automatic way
!  to  check and deal with the status variable coming out of a procedure or
!  set due to an error in part of a process/calculation.
!
!     It will check the statid_chck and will just return if it is equal to
!  zero (no error). If it is different than zero, it will check if argument
!  statid_comm is present so as to copy there information about the problem
!  (statid_pass) and let another part of the program handle it. Please do
!  notice that when returning to the calling procedure, you will still have
!  to handle the problem there (check the example below for the simplest way
!  to do this).
!
!     Finally, if a problem is detected and no statid_comm is provided, it
!  will abort the program with an error message with information about the
!  last place that couldn't deal with the error.
!
!------------------------------------------------------------------------------!
!  ARGUMENTS DESCRIPTION:
!
!  caller_name: Should contain the name of the subroutine calling this check
!               so as to help identify where the error is ocurring.
!
!  statid_chck: The status variable out of a procedure recently called or
!               explicitly set during some calculations. This is the one
!               to be checked against 0 to detect any problem.
!               
!  statid_pass: A numerical value that identifies the location inside the
!               calling subroutine where the error is ocurring. This is the
!               "status" that will be propagated upstream to statid_comm, so
!               that the problem can be dealt with elsewhere (so it should
!               serve to identify what the problem is about).
!
!  statid_comm: If this is present, the content of statid_pass will be then
!               copied here so that another part of the program may deal
!               with this issue. If not, then an error message with some
!               information will be printed and the program will stop.
!
!------------------------------------------------------------------------------!
!  USAGE EXAMPLE:
!
!    subroutine example( arg_A1, ... , arg_AN, extern_stat )
!==>    use liosubs, only: catch_error, [...]
!       implicit none
!       [...]
!       integer, intent(out), optional :: extern_stat
!       integer                        :: intern_stat
!       [...]
!       call subtest( arg_B1, ... , arg_BN, intern_stat )
!==>    catch_error( "example", intern_stat, 11, extern_stat )
!==>    if ( intern_stat /= 0 ) return
!       [...]
!    end subroutine
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   character(len=*), intent(in)            :: caller_name
   integer,          intent(in)            :: statid_chck
   integer,          intent(in)            :: statid_pass
   integer,          intent(out), optional :: statid_comm

   if ( statid_chck == 0) return !  ELSE...

   if ( present(statid_comm) ) then
      statid_comm = statid_pass

   else
      print "(A)",    ""
      print "(A)",    "CRITICAL ERROR OCURRED"
      print "(2A)",   " * Inside Procedure:  ", adjustl(caller_name)
      print "(A,I8)", " * Error place ID:    ", statid_pass
      print "(A,I8)", " * Error cause ID:    ", statid_chck
      print "(A)",    ""
      print "(A)",    "If you can't find out what caused this error, please"
      print "(A)",    "report it back to us through our github page:"
      print "(A)",    ""
      print "(A)",    "   https://github.com/MALBECC/lio"
      print "(A)",    ""
      print "(A)",    "Aborting run..."
      stop

   end if

end subroutine catch_error
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
