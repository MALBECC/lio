!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine catch_iostat(source_name,action_id,unit_id,iost)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  character(len=*),intent(in) :: source_name
  integer,intent(in)          :: action_id
  integer,intent(in)          :: unit_id
  integer,intent(in)          :: iost
!
  character(len=50) :: action_description
  if (iost.ne.0) then
!
!
!
! Identify the operation that caused the error
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
!
!
!
! Write the message
!--------------------------------------------------------------------!
    print '(A)',    ''
    print '(A)',    &
      '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print '(A)',    'IOSTAT CRITICAL ERROR'
    print '(A,I5)', '  IOstat: ', iost
    print '(A,I5)', '  Unit:   ', unit_id
    print '(2A)',   '  Source: ', adjustl(source_name)
    print '(2A)',   '  Action: ', adjustl(action_description)
    print '(A)',    ''
    print '(A)',    ''
    stop
!
!
!
  endif
  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
