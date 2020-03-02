subroutine shift_diag_ON (this, shift_val, offset)

   implicit none
   class(operator), intent(inout) :: this
   LIODBLE         , intent(in)    :: shift_val
   integer        , optional, intent(in) :: offset
   integer :: icount, start = 1

   if (present(offset)) start = offset
   do icount = start, size(this%data_ON,1)
      this%data_ON(icount,icount) = this%data_ON(icount,icount) + shift_val
   enddo

end subroutine shift_diag_ON
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine shift_diag_AO (this, shift_val, offset)

   implicit none
   class(operator), intent(inout) :: this
   LIODBLE         , intent(in)    :: shift_val
   integer        , optional, intent(in) :: offset
   integer :: icount, start = 1

   if (present(offset)) start = offset
   do icount = start, size(this%data_AO,1)
      this%data_AO(icount,icount) = this%data_AO(icount,icount) + shift_val
   enddo

end subroutine shift_diag_AO