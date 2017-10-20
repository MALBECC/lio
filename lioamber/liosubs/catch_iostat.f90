!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine catch_iostat(iost,unit_id,caller_name,action_id)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in)          :: iost
  integer,intent(in)          :: unit_id
  character(len=*),intent(in) :: caller_name
  integer,intent(in)          :: action_id
  character(len=50)           :: action_description


! If there is no error, return to the caller
!--------------------------------------------------------------------!
  if (iost.eq.0) return


! First translate to words the kind of activity being performed
!--------------------------------------------------------------------!
  select case (action_id)
    case (1)
      action_description='Generic inquire operation'
    case (2)
      action_description='Generic open operation'
    case (3)
      action_description='Generic close operation'
    case (4)
      action_description='Generic backspace/rewind operation'
    case (5)
      action_description='Generic read operation'
    case (6)
      action_description='Generic write operation'
    case default
      action_description='Activity unknown...'
  endselect


! Write the error message and stop the program
!--------------------------------------------------------------------!
  print '(A)',    ''
  print '(A)',    &
  '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print '(A)',    'IOSTAT CRITICAL ERROR'
  print '(A,I5)', '  IOstat: ', iost
  print '(A,I5)', '  Unit:   ', unit_id
  print '(2A)',   '  Source: ', adjustl(caller_name)
  print '(2A)',   '  Action: ', adjustl(action_description)
  print '(A)',    ''
  print '(A)',    ''
  stop; end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
