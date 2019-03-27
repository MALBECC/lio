!##############################################################################!
module keywords
   implicit none
   character(len=80) :: output_type    = "full_ground"
   character(len=80) :: output_xyz     = "lioinp_coords"
   character(len=80) :: output_rst     = "lioinp_restart"
   character(len=80) :: inpname_fchk   = "input.fchk"

   namelist /inputgen/ output_type &
   &, output_xyz, output_rst &
   &, inpname_fchk

contains

!##############################################################################!

subroutine reads_keywords( file_name )
   implicit none
   character(len=80), intent(in) :: file_name
   integer                       :: int_stat

   if ( len_trim(file_name) == 0 ) then
      print*, "WARNING IN: read_keywords"
      print*, "No input file provided, using defaults."
      print*; return
   end if

   open( unit= 200, file= file_name, iostat= int_stat )
   if ( int_stat /= 0 ) then
      print*, "WARNING IN: read_keywords"
      print*, "Problem when opening keywords file..."
      print*, "iostat = ", int_stat
      print*; return
   end if

   read( unit= 200, nml= inputgen, iostat= int_stat )
   if ( int_stat > 0 ) then
      print*, "WARNING IN: read_keywords"
      print*, "Problem when reading keywords, using defaults."
      print*, "iostat = ", int_stat
      print*; return
   else if ( int_stat < 0 ) then
      print*, "WARNING IN: read_keywords"
      print*, "Input keywords not found, using defaults."
      print*, "iostat = ", int_stat
      print*; return
   end if

end subroutine reads_keywords

!##############################################################################!

subroutine write_keywords( file_unit )
   implicit none
   integer, intent(in) :: file_unit
   integer :: int_stat

   write( unit=file_unit, fmt='(A)' ) &
   & "------------------------------------------------------------"
   write( unit=file_unit, nml=inputgen, iostat=int_stat )
   if ( int_stat /= 0 ) then
      print*, "ERROR IN: write_keywords"
      print*, "Problem when writing keywords..."
      print*, "iostat = ", int_stat
      stop
   end if
   write( unit=file_unit, fmt='(A)' ) &
   & "------------------------------------------------------------"
   write( unit=file_unit, fmt='(A)' ) ""

end subroutine write_keywords

end module keywords
!##############################################################################!
