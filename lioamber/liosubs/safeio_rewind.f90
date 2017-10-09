!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine safeio_rewind( file_unit, return_stat )
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   integer, intent(in)            :: file_unit
   integer, intent(out), optional :: return_stat

   logical                        :: unit_opened
   integer                        :: mystat
   character(len=20), parameter   :: myname="safeio_rewind"


!  Sanity Check
   mystat = 0
   inquire( unit = file_unit, opened = unit_opened, iostat = mystat )
   call catch_error( myname, mystat, 1, return_stat )
   if ( mystat /= 0 ) return

   if (.not.unit_opened) mystat=1
   call catch_error( myname, 1, 2, return_stat )
   if ( mystat /= 0 ) return


!  Actual Rewind
   mystat = 0
   rewind( unit = file_unit, iostat = mystat )
   call catch_error( myname, mystat, 3, return_stat )
   if ( mystat /= 0 ) return


   if ( present(return_stat) ) return_stat = 0
end subroutine safeio_rewind
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
