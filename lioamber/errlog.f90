!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module errlog

   implicit none
   integer                        :: errlog_count = 0
   integer,           allocatable :: errlog_statids(:)
   integer,           allocatable :: errlog_infoids(:)
   character(len=80), allocatable :: errlog_sources(:)

contains
!------------------------------------------------------------------------------!
! GENERAL DESCRIPTION
!
!     This module was designed to keep track of errors and to help identify
!  their cause and exact location in the code. To do so, the following arrays
!  are used:
!
!     errlog_statids: List of the status that describe the cause of all the
!                     errors that occurred.
!
!     errlog_infoids: List of numbers that help identify the context and/or
!                     location of the error.
!
!     errlog_sources  => List of the subroutines where the errors ocurred.
!
!
!     Two subroutines are provided to control this module, and there are two
!  ways of using it. If you have optional status arguments in your procedures,
!  you can call errlog_Checkin with those optional arguments and whenever a
!  problem araises, the first subroutine that was not provided such an
!  argument will unload the whole error chain:
!
!
!     subroutine example_sub( argA1, ..., argAN, extern_stat )
!        (...)
!        use errlog, only: errlog_Checkin
!        (...)
!        integer, intent(out), optional :: intern_stat
!        (...)
!        call internal_sub( argB1, ..., argBN, intern_stat )
!   =>   call errlog_Checkin("example_sub", intern_stat, [IDNUM], extern_stat)
!   =>   if ( intern_stat /= 0 ) return
!        (...)
!     end subroutine
!
!
!     The second posibility, if you don't have optional status argument in
!  one of your subroutines, is to call errlog_Checkin with a dummy variable
!  in the position of extern_stat
!
!     Notice how errlog_Checkin has the same argument structure as the more
!  simple catch_error (see module liosubs). This allows programers to easily
!  switch from the simple one to this one by using the following alias:
!
!     use liosubs, only: catch_error                     ! (remove this)
!     use errlog,  only: errlog_Checkin => catch_error   ! (add this)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine errlog_Checkin( caller_name, statid, infoid, passid )

   implicit none
   character(len=*), intent(in)            :: caller_name
   integer,          intent(in)            :: statid
   integer,          intent(in)            :: infoid
   integer,          intent(out), optional :: passid

   integer                        :: nn, pi, pf
   integer,           allocatable :: tmplog_statids(:)
   integer,           allocatable :: tmplog_infoids(:)
   character(len=80), allocatable :: tmplog_sources(:)


   if ( present(passid) ) passid=0
   if ( statid == 0 ) return


   if ( errlog_count > 0 ) then
      allocate( tmplog_statids(errlog_count) )
      allocate( tmplog_infoids(errlog_count) )
      allocate( tmplog_sources(errlog_count) )

      tmplog_statids = errlog_statids
      tmplog_infoids = errlog_infoids
      tmplog_sources = errlog_sources

      deallocate( errlog_statids )
      deallocate( errlog_infoids )
      deallocate( errlog_sources )
   end if


   errlog_count = errlog_count + 1
   allocate( errlog_statids(errlog_count) )
   allocate( errlog_infoids(errlog_count) )
   allocate( errlog_sources(errlog_count) )


   pi = 1
   pf = min( 80, len(caller_name) )
   errlog_sources(1) = caller_name( pi : pf )
   errlog_infoids(1) = infoid
   errlog_statids(1) = statid
   do nn = 2, errlog_count
      errlog_sources(nn) = tmplog_sources(nn-1)
      errlog_infoids(nn) = tmplog_infoids(nn-1)
      errlog_statids(nn) = tmplog_statids(nn-1)
   end do


   if ( present(passid) ) then
      passid = infoid

   else
      print "(A)", ""
      print "(A)", "CRITICAL ERROR DETECTED!"
      call errlog_Printlog()
      print "(A)", "Aborting run..."
      print "(A)", ""
      stop

   end if
end subroutine errlog_Checkin



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine errlog_Checkout( last_errid )
   integer, intent(out), optional :: last_errid


   if ( present(last_errid) ) last_errid=0
   if ( errlog_count == 0 ) return


   if ( present(last_errid) ) then
      last_errid = errlog_infoids(1)

   else
      print "(A)", ""
      print "(A)", "CRITICAL ERROR DETECTED!"
      call errlog_Printlog()
      print "(A)", "Aborting run..."
      print "(A)", ""
      stop

   end if
end subroutine errlog_Checkout



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine errlog_Printlog()
   integer :: nn
   print "(A)", "Printing log from last to first ..."
   print "(A)", ""
   do nn = 1, errlog_count
      print "(A,I8)", "ERROR NUMBER: ", errlog_count-nn
      print "(2A)",   "   * Inside Procedure:  ", adjustl( errlog_sources(nn) )
      print "(A,I8)", "   * Location ID:       ", errlog_infoids(nn)
      print "(A,I8)", "   * Error Status:      ", errlog_statids(nn)
   end do
   print "(A)", ""
   print "(A)", "If you don't know and can't find out what caused any of this"
   print "(A)", "these errors, please report it back to us through our github"
   print "(A)", "website: https://github.com/MALBECC/lio"
   print "(A)", ""
end subroutine errlog_Printlog



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
