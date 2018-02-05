!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module lionml_subs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   contains
!
!  Procedures to handle the lionml namelist:
!
!     lionml_Reads: Reads the namelist from a provided file unit.
!     lionml_Write: Writes the namelist to a requested file unit.
!     lionml_Check: Performs consistency checks on the namelist keywords.
!
!     The parameters are self explanatory: file_unit is the unit of an
!  already opened file, and return stat is a status that lets the caller
!  know weather the namelist was correctly read/written/checked or if there
!  was a problem during this. For simplification, the opening and closing of
!  the input/output files must be handled externally.
!
!  TODO: Implement a (simple) logger and use it for the error messages.
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine lionml_Reads( file_unit, extern_stat )

   use lionml_data, only: lionml
   implicit none
   integer, intent(in)            :: file_unit
   integer, intent(out), optional :: extern_stat
   integer                        :: intern_stat

   intern_stat = 0
   rewind( unit = file_unit, iostat = intern_stat )
   if ( intern_stat /= 0 ) then
      print*,"Can't rewind lionml file. Using lionml defaults."
      print*,"iostat = ", intern_stat
      if ( present(extern_stat) ) extern_stat = 1
      return
   end if

   intern_stat = 0
   read( unit = file_unit, nml = lionml, iostat = intern_stat )
   if ( intern_stat /= 0 ) then
      print*,"Can't find lionml namelist. Using lionml defaults"
      print*,"iostat = ", intern_stat
      if ( present(extern_stat) ) extern_stat = 2
      return
   end if

   if ( present(extern_stat) ) extern_stat = 0
end subroutine lionml_Reads
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine lionml_Write( file_unit, extern_stat )

   use lionml_data, only: lionml
   implicit none
   integer, intent(in)            :: file_unit
   integer, intent(out), optional :: extern_stat
   integer                        :: intern_stat
   logical                        :: unit_is_open

   intern_stat = 0
   inquire( unit = file_unit, opened = unit_is_open, iostat = intern_stat )
   if ( intern_stat /= 0 ) then
      print*,"Can't check if output for lionml namelist is opened."
      print*,"iostat = ", intern_stat
      if ( present(extern_stat) ) extern_stat = 1
      return
   end if

   if (.not.unit_is_open) then
      print*,"Output unit for lionml namelist is not opened."
      print*,"iostat = ", intern_stat
      if ( present(extern_stat) ) extern_stat = 2
      return
   end if

   intern_stat = 0
   write( unit = file_unit, nml = lionml, iostat = intern_stat )
   if ( intern_stat /= 0 ) then
      print*,"Could not write lionml in output."
      print*,"iostat = ", intern_stat
      if ( present(extern_stat) ) extern_stat = 3
      return
   end if

   if ( present(extern_stat) ) extern_stat = 0
end subroutine lionml_Write
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine lionml_Check( extern_stat )
 
   implicit none
   integer, intent(out), optional :: extern_stat
   integer                        :: intern_stat

   intern_stat = 0
   !  PERFORM CHECKS

   if ( present(extern_stat) ) extern_stat = 0
end subroutine lionml_Check
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
