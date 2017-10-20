!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine read_list_n(file_name,list)
       implicit none
       character(len=*)     :: file_name
       integer,intent(out)  :: list(:)

       integer :: kk,errorid
!------------------------------------------------------------------------------!

       open(unit=1001,file=file_name)

       do kk=1,size(list)
         read(unit=1001,fmt=*,iostat=errorid) list(kk)
       enddo

       if (errorid.ne.0) then
         print*,'Error while reading list of integers.'
         stop
       end if

!------------------------------------------------------------------------------!
       return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
