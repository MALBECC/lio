!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine get_unit(fileunit,fileunit0)
!--------------------------------------------------------------------!
       implicit none
       integer,intent(out) :: fileunit
       integer,intent(in)  :: fileunit0
       integer             :: limitunit,errorid
       logical             :: notfound

       fileunit=fileunit0
       limitunit=999
!      PRESENT ONLY WORKS IN MODULES! FUCK FORTRAN
!       if (present(fileunit0)) fileunit=fileunit0

       notfound=.true.
       do while (notfound)
         fileunit=fileunit+1
         inquire(unit=fileunit,opened=notfound,iostat=errorid)
         if (fileunit.gt.999) notfound=.false.
         if (errorid.ne.0) notfound=.false.
       enddo

       if (fileunit.gt.limitunit) then
         print*,'ALL UNITS OCCUPIED';stop
       endif

       if (errorid.ne.0) then
         print*,'AN ERROR OCCURRED ON UNIT ',fileunit
     >         ,' (iostat=',errorid,')';stop
       endif


       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
