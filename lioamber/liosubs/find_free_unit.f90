!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine find_free_unit(search_start,search_stop,file_unit)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in)  :: search_start
  integer,intent(in)  :: search_stop
  integer,intent(out) :: file_unit
  integer             :: iost
  logical             :: is_open

  iost=0
  is_open=.true.
  file_unit=max(-1,search_start-1)

  do while ((is_open).and.(file_unit.lt.search_stop))
    file_unit=file_unit+1
    inquire(unit=file_unit,opened=is_open,iostat=iost)
    call catch_iostat('find_free_unit',1,file_unit,iost)
  enddo

  if (is_open) then
    print '(A)', ''
    print '(A)', &
    'Error in find_free_unit: Either no free unit available in'
    print '(A)', &
    'range, or range selected is invalid.'
    print '(A,I6)', '  From unit: ',search_start
    print '(A,I6)', '  To unit:   ',search_stop
    print '(A)',    ''
    stop
  endif

  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
