!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Print_matrix_d( nsize1, nsize2, matrix, fname )
   implicit none
   integer         , intent(in) :: nsize1, nsize2
   real*8          , intent(in) :: matrix( nsize1, nsize2 )
   character(len=*), intent(in) :: fname

   integer                     :: ii, jj
   integer                     :: funit
   integer                     :: int_stat
   logical         , save      :: first_call = .true.
   character(len=*), parameter :: strfmt = "(2x,I6,2x,I6,2x,E15.8)"

#  ifndef DEBUGGING
      return
#  endif

   funit = 1234
   if (first_call) then
      first_call = .false.
      open( unit=funit, file=fname, iostat=int_stat )
   else
      open( unit=funit, file=fname, position="append", iostat=int_stat )
   endif


   if ( int_stat /= 0 ) then
      print*, "ERROR: there is something wrong with output file..."
      print*, "* file = ", fname
      print*, "* unit = ", funit
      print*, "* iost = ", int_stat
      print*, ""

   else
      do jj = 1, nsize2
!      do ii = 1, nsize1
         ii=jj
         write( unit=funit, fmt=strfmt ) ii, jj, matrix(ii,jj)
         if ( int_stat /= 0 ) then
            print*, "ERROR: while trying to print matrix element..."
            print*, "* element: ( ", ii, " , ", jj," )"
            print*, "* file:   ", fname
            print*, "* unit:   ", funit
            print*, "* iost =  ", int_stat
            return
         endif
!      enddo
      enddo

   endif

   if ( int_stat == 0 ) write( unit=funit, fmt=* ) ""
   close( unit=funit )

end subroutine Print_matrix_d
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
