!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  include "../find_free_unit.f90"
  include "../safely_open.f90"
  program safely_open_ut
  implicit none
  integer           :: narg
  character(len=10) :: carg
  integer           :: iost

  integer           :: test_number
  integer           :: unit_io
!
!
!
! Determine chosen test
!------------------------------------------------------------------------------!
  narg=command_argument_count()
  if (narg.lt.1) then
    print*,'Call must contain the number of the test required'
    stop
  endif

  call get_command_argument(1,carg)
  read(unit=carg,fmt=*,iostat=iost) test_number
  if (narg.lt.1) then
    print*,'Could not read the test number from the input argument'
    stop
  endif

  select case (test_number)
!
!
!
! Test 1: Open new file in a given unit
!------------------------------------------------------------------------------!
  case (1)
    print*,''
    print*,'------------------------------------------------------------'
    unit_io=15
    call safely_open(.false.,.true.,'Output-test',unit_io)
    print*,'Test 1 - This should open a new file in unit 15: ',unit_io
    close(unit=unit_io)
    print*,''
!
!
!
! Test 2: Open new file in any unit
!------------------------------------------------------------------------------!
  case (2)
    print*,''
    print*,'------------------------------------------------------------'
    unit_io=0
    call safely_open(.false.,.true.,'Output-test',unit_io)
    print*,'Test 2 - This should open a new file in unit not 0: ',unit_io
    close(unit=unit_io)
    print*,''
!
!
!
! Test 3: Bad request in unit
!------------------------------------------------------------------------------!
  case (3)
    print*,''
    print*,'------------------------------------------------------------'
    unit_io=3
    print*,'Test 3 - This should crash because of bad unit input: '
    call safely_open(.false.,.true.,'Output-test',unit_io)
    close(unit=unit_io)
    print*,''
!
!
!
! Test 4: Open file, then append
!------------------------------------------------------------------------------!
  case (4)
    print*,''
    print*,'------------------------------------------------------------'

    unit_io=0
    call safely_open(.false.,.true.,'Output-test',unit_io)
    print*,'Test 4 - This should open a new file in unit not 0: ',unit_io
    print*,'Test 4 - Now it should write garbage1 into file.'
    write(unit=unit_io,fmt=*) 'garbage1'
    close(unit=unit_io)

    unit_io=0
    call safely_open(.true.,.true.,'Output-test',unit_io)
    print*,'Test 4 - This should open a new file in unit not 0: ',unit_io
    print*,'Test 4 - Now it should overwrite garbage2 into file.'
    write(unit=unit_io,fmt=*) 'garbage2'
    close(unit=unit_io)

    unit_io=100
    call safely_open(.true.,.false.,'Output-test',unit_io)
    print*,'Test 4 - This should open a new file in unit 100: ',unit_io
    print*,'Test 4 - Now it should append garbage3 into file.'
    write(unit=unit_io,fmt=*) 'garbage3'
    close(unit=unit_io)

    print*,''
!
!
!
! Zero prints an end of test message and the default informs that test_number
! is not acceptable
!------------------------------------------------------------------------------!
  case (0)
    print*,''
    print*,'------------------------------------------------------------'
    print*,'End of test'
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
