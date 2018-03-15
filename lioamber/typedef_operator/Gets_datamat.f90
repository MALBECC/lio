subroutine Gets_data_AO (this, Dmat)

   implicit none
   class(operator), intent(in) :: this
   real*8, intent(out)         :: Dmat(:,:)

   Dmat=this%data_AO

end subroutine Gets_data_AO

subroutine Gets_data_ON (this, Dmat)

   implicit none
   class(operator), intent(in) :: this
   real*8, intent(out)         :: Dmat(:,:)

   Dmat=this%data_ON

end subroutine Gets_data_ON
