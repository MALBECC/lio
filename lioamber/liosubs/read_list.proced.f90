!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!subroutine read_list_X( file_name, list )
!   implicit none
!   character(len=*), intent(in)  :: file_name
!   XXXXXXXXXXXXXXXX, intent(out) :: list(:)
!------------------------------------------------------------------------------!
   integer :: kk
   integer :: mystat

   open( unit=1001, file=file_name, iostat=mystat )
   if ( mystat /= 0 ) then
      print*,'CRITICAL ERROR: inside read_list'
      print*,'  Problem when opening unit 1001 to read: ', file_name
      print*,'  IOSTAT = ', mystat
      print*; stop
   endif


   do kk = 1, size(list)

      read( unit=1001, fmt=*, iostat=mystat ) list(kk)
      if ( mystat /= 0 ) then
         print*,'CRITICAL ERROR: inside read_list'
         print*,'  Problem when reading unit 1001, file: ', file_name
         print*,'  List item number ',kk
         print*,'  IOSTAT = ', mystat
         print*; stop
      endif

   enddo
   
!end subroutine read_list_X
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
