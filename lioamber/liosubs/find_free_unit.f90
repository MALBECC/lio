!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine find_free_unit(search_start,search_stop,file_unit)
!
! Finds units that are not opened to open new files.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in)  :: search_start
  integer,intent(in)  :: search_stop
  integer,intent(out) :: file_unit
  integer             :: iost
  logical             :: keep_looking, is_open

  iost=0
  is_open=.true.
  file_unit=max(-1,search_start-1)
  keep_looking=(search_start.le.search_stop)

  do while (keep_looking)
    file_unit=file_unit+1
    inquire(unit=file_unit,opened=is_open,iostat=iost)
    keep_looking=&
    (is_open).and.(iost.eq.0).and.(file_unit.lt.search_stop)
  enddo

  if (iost.ne.0) then
    print*,''
    print*,'Error in find_free_unit:'
    print*,'An error occurred while inquiring unit: ',file_unit
    print*,'(iostat=',iost,')'
    stop
  endif

  if ((file_unit.ge.search_stop).and.(is_open))then
    print*,''
    print*,'Error in find_free_unit:'
    print*,'All units in range were occupied or range is invalid'
    print*,'(from ',search_start,' to ',search_stop,')'
    stop
  endif

  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
