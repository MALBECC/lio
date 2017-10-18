!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!subroutine check_vecsize_n( nsize, vector, vecname, subname )
!   implicit none
!   integer         , intent(in) :: nsize
!   XXXXXXXXXXXXXXXX, intent(in) :: vector(:)
!   character(len=*), intent(in) :: vecname
!   character(len=*), intent(in) :: subname
!------------------------------------------------------------------------------!
   if ( size(vector) /= nsize ) then
      print*,'FATAL ERROR: Size problem incompativility!'
      print*,'  Inside Procedure: ', subname
      print*,'  Vector Name:      ', vecname
      print*,'  Expected size:    ', nsize
      print*,'  Actual size:      ', size(vector)
      print*; stop
   end if
!------------------------------------------------------------------------------!
!end subroutine check_vecsize_X
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
