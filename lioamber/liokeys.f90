!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module liokeys
!
!    This module contains a set of keywords/options that the user can provide
! in order to customize the calculations. A description of each with posible
! values should be included within each definition. Notice how the module and
! this file are called "lio_keywords", whereas the actual namelist is called
! liokeys,
!
!    Some procedures are provided as well: checknml_liokeys allows to do a
! consistency check for all related variables, whereas the readnml_liokeys
! and writenml_liokeys allows the whole list to be read from or writen to a
! given file or unit (two versions exist of each):
!
!    checknml_liokeys( return_stat )
!    readnml_liokeys( file_unit, return_stat )
!    readnml_liokeys( file_name, return_stat )
!    writenml_liokeys( file_unit, return_stat )
!    writenml_liokeys( file_name, return_stat )
!
!    The parameters are self explanatory: file_unit is the unit of an already
! opened file, file_name is the name of a non-opened file, you can chose which
! one to provide. The parameter return_stat is optional and provides a way to
! handle any error (if not provided, errors will halt the program).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none

   integer :: runtype = 0  ! for SCF 
!                     = 1  ! for TD
!                     = 2  ! for Ehrenfest returning forces = 0
!                     = 3  ! for Ehrenfest returning correct forces
!
!
!  External Electrical Field
!------------------------------------------------------------------------------!
   logical :: extEfld_is_on     = .false.
   logical :: extEfld_is_light  = .false.
   logical :: extEfld_is_finite = .false.

   real*8  :: extEfld_wavelength = 0.0d0 ! in nm
   real*8  :: extEfld_timeamp    = 0.0d0 ! in ps
   real*8  :: extEfld_timepos    = 0.0d0 ! in ps

   real*8  :: extEfld_ampx = 0.0d0 ! in au
   real*8  :: extEfld_ampy = 0.0d0 ! in au
   real*8  :: extEfld_ampz = 0.0d0 ! in au
!
!
!  Namelist definition and interfaces
!------------------------------------------------------------------------------!
   namelist /liokeys/ &
   &  runtype                                                                   &
   &, extEfld_is_on, extEfld_is_light, extEfld_is_finite                        &
   &, extEfld_wavelength, extEfld_timeamp, extEfld_timepos                      &
   &, extEfld_ampx, extEfld_ampy, extEfld_ampz

   interface readnml_liokeys
      module procedure readnml_liokeys_fn
      module procedure readnml_liokeys_fu
   end interface readnml_liokeys

   interface writenml_liokeys
      module procedure writenml_liokeys_fn
      module procedure writenml_liokeys_fu
   end interface writenml_liokeys

contains
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine checknml_liokeys( return_stat )

   use liosubs, only: catch_error
   implicit none
   integer, intent(out), optional :: return_stat

   integer                        :: mystat
   character(len=*), parameter    :: myname="checknml_liokeys"

   mystat = 0
!  PERFORM CHECKS
   call catch_error( myname, mystat, mystat, return_stat )

   if ( present(return_stat) ) return_stat = 0
end subroutine checknml_liokeys
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine readnml_liokeys_fu( file_unit, return_stat )

   use liosubs, only: catch_error, safeio_rewind
   implicit none
   integer, intent(in)            :: file_unit
   integer, intent(out), optional :: return_stat

   integer                        :: mystat
   character(len=*), parameter    :: myname="readnml_liokeys_fu"
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
   read( unit = file_unit, nml = liokeys, iostat = mystat )
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
end subroutine readnml_liokeys_fu
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine readnml_liokeys_fn( file_name, return_stat )

   use liosubs, only: catch_error, safeio_open
   implicit none
   character(len=*), intent(in)   :: file_name
   integer, intent(out), optional :: return_stat

   integer                        :: file_unit
   integer                        :: mystat
   character(len=*), parameter    :: myname="readnml_liokeys_fn"
!
!
!  Open the file in an available unit
!------------------------------------------------------------------------------!
   mystat = 0
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
   call readnml_liokeys_fu( file_unit, mystat )
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
end subroutine readnml_liokeys_fn
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine writenml_liokeys_fu( file_unit, return_stat )

   use liosubs, only: catch_error
   implicit none
   integer, intent(in)            :: file_unit
   integer, intent(out), optional :: return_stat

   logical                        :: unit_opened
   integer                        :: mystat
   character(len=*), parameter    :: myname="writenml_liokeys_fu"
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
   write( unit = file_unit, nml = liokeys, iostat = mystat )
   call catch_error( myname, mystat, 3, return_stat )
   if ( mystat /= 0 ) return


   if ( present(return_stat) ) return_stat = 0
end subroutine writenml_liokeys_fu
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine writenml_liokeys_fn( file_name, return_stat )

   use liosubs, only: catch_error, safeio_open
   implicit none
   character(len=*), intent(in)   :: file_name
   integer, intent(out), optional :: return_stat

   integer                        :: file_unit
   integer                        :: mystat
   character(len=*), parameter    :: myname="writenml_liokeys_fu"
!
!
!  Open, write, close, pretty clear which is which
!------------------------------------------------------------------------------!
   mystat = 0
   call safeio_open( file_unit, file_name, 3, mystat )
   call catch_error( myname, mystat, 1, return_stat )
   if ( mystat /= 0 ) return

   mystat = 0
   call writenml_liokeys_fu( file_unit, mystat )
   call catch_error( myname, mystat, 2, return_stat )
   if ( mystat /= 0 ) return

   mystat = 0
   close( unit = file_unit, iostat = mystat )
   call catch_error( myname, mystat, 3, return_stat )
   if ( mystat /= 0 ) return

   if ( present(return_stat) ) return_stat = 0
end subroutine writenml_liokeys_fn
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
