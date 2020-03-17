!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Gets_smat( this, Smat )

   implicit none
   class(sop), intent(in)  :: this
   LIODBLE    , intent(out) :: Smat( this%Nbasis, this%Nbasis )

   Smat = this%Smat

end subroutine Gets_smat
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
