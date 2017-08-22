!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine safeio_open( file_unit, file_name, open_mode, return_stat )
!
!  open_mode = ...
!  >
!  > = +-1   Opens existing file to read. 
!  >
!  > = +-2   Opens/creates file to write (appending existing).
!  >
!  > = +-3   Opens/creates file to write (overwritting existing).
!  >
!  > (positive numbers for formatted, negative for binary)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   integer,          intent(inout)          :: file_unit
   character(len=*), intent(in)             :: file_name
   integer,          intent(in)             :: open_mode
   integer,          intent(out),  optional :: return_stat

   logical                      :: file_opened
   character(len=20)            :: file_format
   integer                      :: mystat
   character(len=20), parameter :: myname="safeio_open"
!
!
!
! (1) Initialization and input mode check
!------------------------------------------------------------------------------!
   if ( present(return_stat) ) return_stat = 0
   mystat = 0

   if ( open_mode > 0 ) then
      file_format = "formatted"
   else if ( open_mode > 0 ) then
      file_format = "unformatted"
   else
      mystat = 1
   endif

   if ( abs(open_mode) > 3 ) then
      mystat = open_mode
   endif

   call catch_error( myname, mystat, 1, return_stat )
   if ( mystat /= 0 ) return
!
!
!
! (2) Check that the file is not opened and there is available units
!------------------------------------------------------------------------------!
   mystat = 0
   inquire( file = file_name, opened = file_opened, iostat = mystat )
   call catch_error( myname, mystat, 2, return_stat )
   if ( mystat /= 0 ) return

   if (file_opened) mystat = 1
   call catch_error( myname, mystat, 3, return_stat )
   if ( mystat /= 0 ) return

   if ( file_unit < 10 ) file_unit = 10
   call find_free_unit( file_unit, 1000, mystat )
   call catch_error( myname, mystat, 4, return_stat )
   if ( mystat /= 0 ) return
!
!
!
! (3) Choose appropriate way of opening and do so
!------------------------------------------------------------------------------!
   mystat = 0
   select case (open_mode)

      case (-1,1)
         open( unit = file_unit,   file = file_name,    iostat = mystat,       &
             & form = file_format, status = "old",      action = "read" )

      case (-2,2)
         open( unit = file_unit,   file = file_name,    iostat = mystat,       &
             & form = file_format, position = "append", action = "write" )

      case (-3,3)
         open( unit = file_unit,   file = file_name,    iostat = mystat,       &
             & form = file_format, status = "replace",  action = "write" )

   end select
   call catch_error( myname, mystat, 5, return_stat )
   if ( mystat /= 0 ) return
!
!
!
!------------------------------------------------------------------------------!
end subroutine safeio_open
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
