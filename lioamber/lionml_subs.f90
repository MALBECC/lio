!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%  LIONML_SUBS.F90  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains procedures to handle both lio and lionml namelists. It    !
! includes the following subroutines:                                          !
! * lionml_read: Reads both lio and lionml namelists from input file.          !
! * lionml_write: Prints both lio and lionml namelists to standard output.     !
! * lionml_check: Performs consistency checks on the namelist keywords. (TO-DO)!
!                                                                              !
! In addition, the following subroutines are meant only accessed internally:   !
! * lionml_write_dull: Prints namelists in a straightforward manner.           !
! * lionml_write_style: Prints namelists in a fancy manner.                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module lionml_subs
   implicit none
contains

!  The parameters are self explanatory: file_unit is the unit of an already
!  opened file, and return stat is a status that lets the caller know wether
!  the namelist was correctly read/written/checked or if there was a problem
!  during these processes. For sake of simplification, the opening and closing
!  of input/output files must be handled externally.
!  TODO: Implement a (simple) logger and use it for the error messages.
subroutine lionml_check(extern_stat)
   implicit none
   integer, intent(out), optional :: extern_stat
   integer                        :: intern_stat

   intern_stat = 0
   if ( present(extern_stat) ) extern_stat = 0
   return
end subroutine lionml_check

subroutine lionml_read(file_unit, extern_stat )

   use lionml_data, only: lionml, lio
   use fileio_data, only: verbose
   implicit none
   integer, intent(in)            :: file_unit
   integer, intent(out), optional :: extern_stat
   integer                        :: intern_stat

   ! Old lio namelist.
   intern_stat = 0
   rewind( unit = file_unit, iostat = intern_stat )
   if ( intern_stat /= 0 ) then
      write(*,'(A)') &
         "Cannot rewind LIO input file. Using defaults for namelist lio."
      if (verbose .gt. 3) write(*,'(A,I4)') "iostat = ", intern_stat
      if ( present(extern_stat) ) extern_stat = 1
      return
   end if

   intern_stat = 0
   read( unit = file_unit, nml = lio, iostat = intern_stat )
   if ( intern_stat /= 0 ) then
      write(*,'(A)') &
         "Cannot find lio namelist. Using defaults for namelist lio."
      if (verbose .gt. 3) write(*,'(A,I4)') "iostat = ", intern_stat
      if ( present(extern_stat) ) extern_stat = 2
      return
   end if

   ! New lionml namelist.
   intern_stat = 0
   rewind( unit = file_unit, iostat = intern_stat )
   if ( intern_stat /= 0 ) then
      write(*,'(A)') &
         "Cannot rewind LIO input file. Using defaults for namelist lionml."
      if (verbose .gt. 3) write(*,'(A,I4)') "iostat = ", intern_stat
      if ( present(extern_stat) ) extern_stat = 1
      return
   end if

   intern_stat = 0
   read( unit = file_unit, nml = lionml, iostat = intern_stat )
   if ( intern_stat /= 0 ) then
      write(*,'(A)') &
         "Cannot find lionml namelist. Using defaults for namelist lionml."
      if (verbose .gt. 3) write(*,'(A,I4)') "iostat = ", intern_stat
      if ( present(extern_stat) ) extern_stat = 2
      return
   end if

   if ( present(extern_stat) ) extern_stat = 0
   return
end subroutine lionml_read


subroutine lionml_write(extern_stat)
   use fileio_data, only: get_style
   implicit none
   integer, intent(out), optional :: extern_stat
   logical                        :: my_style

   call get_style(my_style)
   if (my_style) then
      call lionml_write_style()
   else
      call lionml_write_dull()
   endif

   return
end subroutine lionml_write

subroutine lionml_write_dull()


end subroutine lionml_write_dull



subroutine lionml_write_style()


end subroutine lionml_write_style

end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
