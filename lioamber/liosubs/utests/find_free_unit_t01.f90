!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  include "../find_free_unit.f90"
  program find_free_unit_t01
  implicit none
  integer           :: test_number
  integer           :: unit_found
  integer           :: iost
  integer           :: narg
  character(len=10) :: carg
!
!
! Determine chosen test
!------------------------------------------------------------------------------!
  narg=command_argument_count()
  if (narg.lt.1) then
    print*,'Error in find_free_unit_t01'
    print*,'Call must contain the number of the test required'
    stop
  endif
  call get_command_argument(1,carg)

  read(unit=carg,fmt=*,iostat=iost) test_number
  if (narg.lt.1) then
    print*,'Error in find_free_unit_t01'
    print*,'Could not read the test number from the input argument'
    stop
  endif
  select case (test_number)
!
!
! Test 1: Chooses first option
!------------------------------------------------------------------------------!
  case (1)
    print*,''
    print*,'------------------------------------------------------------'
    call find_free_unit(20,22,unit_found)
    print*,'Test 1 - This number should be 20: ',unit_found
    print*,''
!
!
! Test 2: Skips opened file
!------------------------------------------------------------------------------!
  case (2)
    print*,''
    print*,'------------------------------------------------------------'
    open(file='Output-empty20.o',unit=20)
    call find_free_unit(20,22,unit_found)
    print*,'Test 2 - This number should be 21: ',unit_found
    close(unit=20)
    print*,''
!
!
! Test 3: Error message if all files are open
!------------------------------------------------------------------------------!
  case (3)
    print*,''
    print*,'------------------------------------------------------------'
    open(file='Output-empty20.o',unit=20)
    open(file='Output-empty21.o',unit=21)
    open(file='Output-empty22.o',unit=22)
    print*,'Test 3 - Error because units 20 to 22 are occuppied:...'
    call find_free_unit(20,22,unit_found)
    close(unit=20)
    close(unit=21)
    close(unit=22)
    print*,''
!
!
! Test 3: Error message if range doesn't make sense
!------------------------------------------------------------------------------!
  case (4)
    print*,''
    print*,'------------------------------------------------------------'
    print*,'Test 4 - Error because of bad range:...'
    call find_free_unit(22,20,unit_found)
!
!
! Default and end
!------------------------------------------------------------------------------!
  case (0)
    print*,''
    print*,'------------------------------------------------------------'
    print*,'End of testing'
    print*,''
    print*,''

  case default
    print*,''
    print*,'------------------------------------------------------------'
    print*,'Wrongh test id: ',test_number
    print*,''
  endselect

  return; end program
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
