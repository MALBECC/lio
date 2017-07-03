!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine find_free_unit( free_unit, last_unit, return_stat )
!
!     This subroutine automatizes the process of finding a free available
!  unit to open a file into. It receives three arguments:
!
!  free_unit (inout):
!     On input, it contains the first unit that the subroutine will try.
!     On output, it contains the last unit that the subroutine has tried
!     (if no error occurrs, that last unit is one that is free and available
!     to use). This subroutine will never use a unit lower that 10 (in which
!     case, it ignores this argument).
!
!  last_unit (in):
!     The last unit the subroutine will try out.
!
!  resturn_stat (optional out):
!     Contains information about any error that might have occurred during
!     this procedure. If nothing went wrongh, it returns 0. If it is not
!     present, any error will halt the program.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   integer, intent(inout)          :: free_unit
   integer, intent(in)             :: last_unit
   integer, intent(out),  optional :: return_stat

   integer                      :: mystat
   character(len=20), parameter :: myname="find_free_unit"

   integer :: first_try
   integer :: final_try
   integer :: test_unit
   logical :: unit_is_open


   mystat = 0
   first_try = max( free_unit, 10 )
   final_try = last_unit
   if ( first_try > final_try ) mystat = 1
   call catch_error( myname, mystat, 1, return_stat )
   if ( mystat /= 0 ) return

   do test_unit = first_try, final_try
      unit_is_open = .true.
      inquire( unit = test_unit, opened = unit_is_open, iostat = mystat )
      if (.not.unit_is_open) exit
   end do

   if (unit_is_open) then
      mystat = 1
      test_unit = -test_unit
   end if
   call catch_error( myname, mystat, test_unit, return_stat )
   if ( mystat /= 0 ) return

   free_unit = test_unit
   if ( present(return_stat) ) return_stat = 0
end subroutine find_free_unit
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
