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

   logical           :: file_exists
   logical           :: file_opened

   character(len=20) :: file_action
   character(len=20) :: file_format
   character(len=20) :: file_status

   logical                      :: can_create, can_replace
   integer                      :: mystat
   character(len=20), parameter :: myname="safeio_open"



!  (1) Apply setup according to option selected
!------------------------------------------------------------------------------!
   mystat = 0
   select case (open_mode)
      case (1,-1)
         file_action = "read"
         can_create  = .false.
         can_replace = .false.

      case (2,-2)
         file_action = "write"
         can_create  = .true.
         can_replace = .false.

      case (3,-3)
         file_action = "write"
         can_create  = .true.
         can_replace = .true.

      case default
         mystat = open_mode

   end select
   call catch_error( myname, mystat, 1, return_stat )
   if ( mystat /= 0 ) return

   if ( open_mode > 0 ) then
      file_format = "formatted"
   else
      file_format = "unformatted"
   end if


!  (2) Inquire if file exists and if its opened and adapt to it
!------------------------------------------------------------------------------!
   mystat = 0
   inquire( file = file_name,   exist = file_exists, opened = file_opened, &
          & iostat = mystat )
   call catch_error( myname, mystat, 2, return_stat )
   if ( mystat /= 0 ) return


   if (file_opened) then
      call catch_error( myname, 1, 3, return_stat )
      return
   end if


   if (file_exists) then
      if (can_replace) then
         file_status = "replace"
      else
         file_status = "old"
      end if

   else
      file_status = "new"
      if (.not.can_create) then
        call catch_error( myname, mystat, 4, return_stat )
        return
      end if

   end if



!  (3) Find a free unit and open the file there
!------------------------------------------------------------------------------!
   mystat = 0
   call find_free_unit( file_unit, 1000, mystat )
   call catch_error( myname, mystat, 5, return_stat )
   if ( mystat /= 0 ) return


   mystat = 0
   open( unit = file_unit,   file = file_name,     iostat = mystat,          &
       & form = file_format, status = file_status, action = file_action )
   call catch_error( myname, mystat, 6, return_stat )
   if ( mystat /= 0 ) return


   if ( present(return_stat) ) return_stat = 0
end subroutine safeio_open
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
