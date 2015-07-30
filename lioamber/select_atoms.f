!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       pure subroutine select_atoms
     > (N,chosen_group,atom_group,atom_selection)

       implicit none
       integer,intent(in)   :: N,chosen_group
       integer,intent(in)   :: atom_group(N)
       integer,intent(out)  :: atom_selection(N)

       integer :: kk
!------------------------------------------------------------------------------!

       do kk=1,N
         atom_selection(kk)=0
         if (atom_group(kk).EQ.chosen_group) atom_selection(kk)=1
       enddo

!------------------------------------------------------------------------------!
       return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
