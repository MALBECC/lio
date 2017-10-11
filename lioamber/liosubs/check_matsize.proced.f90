!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!subroutine check_matsize_X( nsize1, nsize2, matrix, matname, subname )
!   implicit none
!   integer          , intent(in) :: nsize1
!   integer          , intent(in) :: nsize2
!   XXXXXXXXXXXXXXXXX, intent(in) :: matrix(:)
!   character(len=80), intent(in) :: matname
!   character(len=80), intent(in) :: subname
!------------------------------------------------------------------------------!
   if (( size(matrix,1) /= nsize1 ).or.( size(matrix,2) /= nsize2 )) then
      print*,'FATAL ERROR: Size problem incompativility!'
      print*,'  Inside Procedure: ', subname
      print*,'  Matrix Name:      ', matname
      print*,'  Expected sizes:   ', nsize1, nsize2
      print*,'  Actual sizes:     ', size(matrix,1), size(matrix,2)
      print*; stop
   end if
!------------------------------------------------------------------------------!
!end subroutine check_matsize_X
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
