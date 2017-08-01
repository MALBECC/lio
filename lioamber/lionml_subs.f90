!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module lionml_subs
!
!    Some procedures are provided as well: lionml_Check allows to do a
! consistency check for all related variables, whereas the lionml_Read
! and lionml_Write allows the whole list to be read from or writen to a
! given file or unit (two versions exist of each):
!
!    lionml_Check( return_stat )
!    lionml_Read_fu( file_unit, return_stat )
!    lionml_Read_fn( file_name, return_stat )
!    lionml_Write_fu( file_unit, return_stat )
!    lionml_Write_fn( file_name, return_stat )
!
!    The parameters are self explanatory: file_unit is the unit of an already
! opened file, file_name is the name of a non-opened file, you can chose which
! one to provide. The parameter return_stat is optional and provides a way to
! handle any error (if not provided, errors will halt the program).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none

   interface lionml_Read
      module procedure lionml_Read_fn
      module procedure lionml_Read_fu
   end interface lionml_Read

   interface lionml_Write
      module procedure lionml_Write_fn
      module procedure lionml_Write_fu
   end interface lionml_Write

contains
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine lionml_Check( return_stat )

   use liosubs, only: catch_error
   implicit none
   integer, intent(out), optional :: return_stat

   integer                        :: mystat
   character(len=*), parameter    :: myname="lionml_Check"

   mystat = 0
!  PERFORM CHECKS
   call catch_error( myname, mystat, mystat, return_stat )

   if ( present(return_stat) ) return_stat = 0
end subroutine lionml_Check
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine lionml_Read_fu( file_unit, return_stat )

   use liosubs, only: catch_error, safeio_rewind
   use lionml_data, only: lionml
   implicit none
   integer, intent(in)            :: file_unit
   integer, intent(out), optional :: return_stat

   integer                        :: mystat
   character(len=*), parameter    :: myname="lionml_Read_fu"
!
!
!  Rewind to make sure the namelist has not been already passed
!------------------------------------------------------------------------------!
   mystat = 0
   call safeio_rewind( file_unit, mystat )
   call catch_error( myname, mystat, 1, return_stat )
   if ( mystat /= 0 ) return
!
!
!  Read the namelist (if wasn't found, return error)
!------------------------------------------------------------------------------!
   mystat = 0
   read( unit = file_unit, nml = lionml, iostat = mystat )
   if ( (mystat == 84) .or. (mystat == 85) ) then
!     (case: namelist not found)
      call catch_error( myname, mystat, 2, return_stat )
      return
   else
!     (case: other errors)
      call catch_error( myname, mystat, 3, return_stat )
      if ( mystat /= 0 ) return
   end if


   if ( present(return_stat) ) return_stat = 0
end subroutine lionml_Read_fu
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine lionml_Read_fn( file_name, return_stat )

   use liosubs, only: catch_error, safeio_open
   use lionml_data, only: lionml
   implicit none
   character(len=*), intent(in)   :: file_name
   integer, intent(out), optional :: return_stat

   integer                        :: file_unit
   integer                        :: mystat
   character(len=*), parameter    :: myname="lionml_Read_fn"
!
!
!  Open the file in an available unit
!------------------------------------------------------------------------------!
   mystat = 0
   file_unit = 10
   call safeio_open( file_unit, file_name, 1, mystat )
   if ( mystat == 5 ) then
!     (case: file not found)
      call catch_error( myname, mystat, 1, return_stat )
      return
   else
!     (case: other errors)
      call catch_error( myname, mystat, 2, return_stat )
      if ( mystat /= 0 ) return
   end if
!
!
!  Read the namelist
!------------------------------------------------------------------------------!
   mystat = 0
   call lionml_Read_fu( file_unit, mystat )
   if (mystat == 2) then
!     (case: namelist not found)
      call catch_error( myname, mystat, 3, return_stat )
      return
   else
!     (case: other errors)
      call catch_error( myname, mystat, 4, return_stat )
      if ( mystat /= 0 ) return
   end if
!
!
!  Close the open unit
!------------------------------------------------------------------------------!
   mystat = 0
   close( unit = file_unit, iostat = mystat )
   call catch_error( myname, mystat, 5, return_stat )
   if ( mystat /= 0 ) return


   if ( present(return_stat) ) return_stat = 0
end subroutine lionml_Read_fn
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine lionml_Write_fu( file_unit, return_stat )

   use liosubs, only: catch_error
   use lionml_data, only: lionml
   implicit none
   integer, intent(in)            :: file_unit
   integer, intent(out), optional :: return_stat

   logical                        :: unit_opened
   integer                        :: mystat
   character(len=*), parameter    :: myname="lionml_Write_fu"
!
!
!  Sanity Check: confirm that unit is opened
!------------------------------------------------------------------------------!
   mystat = 0
   inquire( unit = file_unit, opened = unit_opened, iostat = mystat )
   if (.not.unit_opened) then
!     (case: unit was not opened)
      call catch_error( myname, 1, 1, return_stat )
      return
   else
!     (case: other errors)
      call catch_error( myname, mystat, 2, return_stat )
      if ( mystat /= 0 ) return
   end if
!
!
!  Write the namelist
!------------------------------------------------------------------------------!
   mystat = 0
   write( unit = file_unit, nml = lionml, iostat = mystat )
   call catch_error( myname, mystat, 3, return_stat )
   if ( mystat /= 0 ) return


   if ( present(return_stat) ) return_stat = 0
end subroutine lionml_Write_fu
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine lionml_Write_fn( file_name, return_stat )

   use liosubs, only: catch_error, safeio_open
   use lionml_data, only: lionml
   implicit none
   character(len=*), intent(in)   :: file_name
   integer, intent(out), optional :: return_stat

   integer                        :: file_unit
   integer                        :: mystat
   character(len=*), parameter    :: myname="lionml_Write_fn"
!
!
!  Open, write, close, pretty clear which is which
!------------------------------------------------------------------------------!
   mystat = 0
   call safeio_open( file_unit, file_name, 3, mystat )
   call catch_error( myname, mystat, 1, return_stat )
   if ( mystat /= 0 ) return

   mystat = 0
   call lionml_Write_fu( file_unit, mystat )
   call catch_error( myname, mystat, 2, return_stat )
   if ( mystat /= 0 ) return

   mystat = 0
   close( unit = file_unit, iostat = mystat )
   call catch_error( myname, mystat, 3, return_stat )
   if ( mystat /= 0 ) return

   if ( present(return_stat) ) return_stat = 0
end subroutine lionml_Write_fn
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
