!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine read_groups(N,file_name,atom_group)
       implicit none
       integer,intent(in)   :: N
       character(len=*)     :: file_name
       integer,intent(out)  :: atom_group(N)

       integer :: kk
!------------------------------------------------------------------------------!

       open(unit=1001,file=file_name)
       do kk=1,N
         read(unit=1001,fmt=*) atom_group(kk)
       enddo

!------------------------------------------------------------------------------!
       return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
