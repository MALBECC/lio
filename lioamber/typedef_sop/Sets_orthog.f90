!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Sets_orthog( this, method_id, maxval_ld )

   implicit none
   class(sop), intent(inout) :: this
   integer   , intent(in)    :: method_id
   real*8    , intent(in)    :: maxval_ld

   if (method_id /= 0) call this%Gets_orthog_4m &
   &  ( method_id, maxval_ld, this%Xmat, this%Ymat, this%Xtrp, this%Ytrp )

end subroutine Sets_orthog
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
