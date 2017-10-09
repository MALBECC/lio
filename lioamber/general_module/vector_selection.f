!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine vector_selection
     > (chosen_id,item_id,selection)

       implicit none
       integer,intent(in)   :: chosen_id
       integer,intent(in)   :: item_id(:)
       integer,intent(out)  :: selection(:)

       integer :: kk,is_selected

!------------------------------------------------------------------------------!

       if (size(item_id).ne.size(selection)) then
         print*,'vector_selection: Size of arrays is not compatible'
         stop
       end if

       do kk=1,size(item_id)
         is_selected=0
         if (item_id(kk).EQ.chosen_id) is_selected=1
         selection(kk)=is_selected
       end do

!------------------------------------------------------------------------------!

       return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
