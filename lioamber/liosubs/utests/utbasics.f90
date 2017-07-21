!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine get_testid( testid )
   integer, intent(out) :: testid
   integer              :: iost
   integer              :: narg
   character(len=10)    :: carg

   narg = command_argument_count()
   if ( narg /= 1 ) then
      print*,"ERROR IN READ_TESTID:"
      print*,"Call must contain ONE requested test id number."
      stop
   end if

   call get_command_argument( 1, carg )
   read( unit=carg, fmt=*, iostat=iost) testid
   if ( iost /= 0 ) then
      print*,"ERROR IN READ_TESTID:"
      print*,'Could not read the test number from the input argument'
      stop
   end if

end subroutine get_testid
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
