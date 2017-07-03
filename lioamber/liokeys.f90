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
! consistency check for all related variables, whereas the liokeys_Readnml
! and liokeys_Writenml allows the whole list to be read from or writen to a
! given file or unit (two versions exist of each):
!
!    checknml_liokeys( return_stat )
!    liokeys_Readnml( file_unit, return_stat )
!    liokeys_Readnml( file_name, return_stat )
!    liokeys_Writenml( file_unit, return_stat )
!    liokeys_Writenml( file_name, return_stat )
!
!    The parameters are self explanatory: file_unit is the unit of an already
! opened file, file_name is the name of a non-opened file, you can chose which
! one to provide. The parameter return_stat is optional and provides a way to
! handle any error (if not provided, errors will halt the program).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
!
   integer :: ndyn_steps = 0  ! Number of total nuclear dynamic steps
   integer :: edyn_steps = 0  ! Number of total electronic dynamic steps PER
                              ! nuclear dynamic step.
!
!  ndyn == 0 & edyn == 0   =>   Single point
!  ndyn == 0 & edyn /= 0   =>   TD electron dynamic
!  ndyn /= 0 & edyn == 0   =>   BO atomistic dynamic
!  ndyn /= 0 & edyn /= 0   =>   Ehrenfest dynamic
!
   logical :: nullify_forces = .false.
!
!
!  Restarting information
!------------------------------------------------------------------------------!
!     If (rst_filei != ""), the program will use the information there to
!  restart the calculations. If (rst_nfreq != 0), the program will save all
!  the necesary information to restart calculations every rstfreq steps in
!  "rst_fileo".
!
   character(len=80) :: rst_filei = ""
   character(len=80) :: rst_fileo = "liorst.out"
   integer           :: rst_nfreq = 0
!
!
!  External Electrical Field
!------------------------------------------------------------------------------!
!     An external field will be active starting from "stepi" (no external
!  field if stepi<0) to "stepf" (always on if stepf < stepi). The field may
!  have two modulations: an oscilating one (if wavelen /= 0) or a gaussian
!  one (if timeamp != 0; will be a full gaussian if field has a stepf or half
!  a gaussian if not).
!
   integer :: eefld_stepi = 0
   integer :: eefld_stepf = 0
   real*8  :: eefld_timeamp = 0.0d0 ! in ps
   real*8  :: eefld_wavelen = 0.0d0 ! in nm

   real*8  :: eefld_ampx = 0.0d0 ! in au
   real*8  :: eefld_ampy = 0.0d0 ! in au
   real*8  :: eefld_ampz = 0.0d0 ! in au
!
!
!  Namelist definition and interfaces
!------------------------------------------------------------------------------!
   namelist /lionml/ &
   &  ndyn_steps, edyn_steps, nullify_forces                                   &
   &, rst_filei, rst_fileo, rst_nfreq                                          &
   &, eefld_stepi, eefld_stepf, eefld_timeamp, eefld_wavelen                   &
   &, eefld_ampx, eefld_ampy, eefld_ampz
!
   interface liokeys_Readnml
      module procedure liokeys_Readnml_fn
      module procedure liokeys_Readnml_fu
   end interface liokeys_Readnml
!
   interface liokeys_Writenml
      module procedure liokeys_Writenml_fn
      module procedure liokeys_Writenml_fu
   end interface liokeys_Writenml
!
contains
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine liokeys_Check( return_stat )

   use liosubs, only: catch_error
   implicit none
   integer, intent(out), optional :: return_stat

   integer                        :: mystat
   character(len=*), parameter    :: myname="liokeys_Check"

   mystat = 0
!  PERFORM CHECKS
   call catch_error( myname, mystat, mystat, return_stat )

   if ( present(return_stat) ) return_stat = 0
end subroutine liokeys_Check
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine liokeys_Readnml_fu( file_unit, return_stat )

   use liosubs, only: catch_error, safeio_rewind
   implicit none
   integer, intent(in)            :: file_unit
   integer, intent(out), optional :: return_stat

   integer                        :: mystat
   character(len=*), parameter    :: myname="liokeys_Readnml_fu"
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
end subroutine liokeys_Readnml_fu
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine liokeys_Readnml_fn( file_name, return_stat )

   use liosubs, only: catch_error, safeio_open
   implicit none
   character(len=*), intent(in)   :: file_name
   integer, intent(out), optional :: return_stat

   integer                        :: file_unit
   integer                        :: mystat
   character(len=*), parameter    :: myname="liokeys_Readnml_fn"
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
   call liokeys_Readnml_fu( file_unit, mystat )
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
end subroutine liokeys_Readnml_fn
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine liokeys_Writenml_fu( file_unit, return_stat )

   use liosubs, only: catch_error
   implicit none
   integer, intent(in)            :: file_unit
   integer, intent(out), optional :: return_stat

   logical                        :: unit_opened
   integer                        :: mystat
   character(len=*), parameter    :: myname="liokeys_Writenml_fu"
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
end subroutine liokeys_Writenml_fu
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine liokeys_Writenml_fn( file_name, return_stat )

   use liosubs, only: catch_error, safeio_open
   implicit none
   character(len=*), intent(in)   :: file_name
   integer, intent(out), optional :: return_stat

   integer                        :: file_unit
   integer                        :: mystat
   character(len=*), parameter    :: myname="liokeys_Writenml_fn"
!
!
!  Open, write, close, pretty clear which is which
!------------------------------------------------------------------------------!
   mystat = 0
   call safeio_open( file_unit, file_name, 3, mystat )
   call catch_error( myname, mystat, 1, return_stat )
   if ( mystat /= 0 ) return

   mystat = 0
   call liokeys_Writenml_fu( file_unit, mystat )
   call catch_error( myname, mystat, 2, return_stat )
   if ( mystat /= 0 ) return

   mystat = 0
   close( unit = file_unit, iostat = mystat )
   call catch_error( myname, mystat, 3, return_stat )
   if ( mystat /= 0 ) return

   if ( present(return_stat) ) return_stat = 0
end subroutine liokeys_Writenml_fn
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
