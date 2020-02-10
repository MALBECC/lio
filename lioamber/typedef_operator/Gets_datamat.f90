!carlos: these subroutines extract matrices from operator.
subroutine Gets_data_AO (this, Dmat)

   implicit none
   class(operator), intent(in) :: this
   real*8, intent(out)         :: Dmat(:,:)

   Dmat=this%data_AO

end subroutine Gets_data_AO
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Gets_data_ON (this, Dmat)

   implicit none
   class(operator), intent(in) :: this
   real*8, intent(out)         :: Dmat(:,:)

   Dmat=this%data_ON

end subroutine Gets_data_ON
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Gets_dataC_AO (this, Dmat)

   implicit none
   class(operator), intent(in)  :: this
#ifdef TD_SIMPLE
   complex*8, intent(out)       :: Dmat(:,:)
#else
    complex*16, intent(out)     :: Dmat(:,:)
#endif

   Dmat=this%dataC_AO

end subroutine Gets_dataC_AO
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Gets_dataC_ON (this, Dmat)

   implicit none
   class(operator), intent(in) :: this
#ifdef TD_SIMPLE
   complex*8, intent(out)      :: Dmat(:,:)
#else
    complex*16, intent(out)    :: Dmat(:,:)
#endif

   Dmat=this%dataC_ON

end subroutine Gets_dataC_ON
