!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine safely_open(should_exist,from_start,file_name,file_unit)
!
! Opens file_name in file_unit. Looks for an open unit if file_unit
! is 0.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  logical,intent(in)          :: should_exist
  logical,intent(in)          :: from_start
  character(len=*),intent(in) :: file_name
  integer,intent(inout)       :: file_unit

  integer :: iost
  logical :: file_exists
!
!
!
! Check if input unit is compatible
!--------------------------------------------------------------------!
  if (file_unit.lt.10) then
    if (file_unit.eq.0) then
      call find_free_unit(10,1000,file_unit)
    else
      print*,'Error in safely_open:'
      print*,'Unit requested should be greater than 10, but is '&
            ,file_unit
      print*,'(or set to 0 if you want to look for any free unit)'
      print*,''
    endif
  endif
!
!
!
! Check if file requested is present
!--------------------------------------------------------------------!
  if (should_exist) then
    inquire(file=file_name,exist=file_exists,iostat=iost)

    if (.not.file_exists) then
      print*,'Error in safely_open:'
      print*,'The file ',file_name,' does not exist (and it should).'
      print*,''
    endif

    if (iost.ne.0) then
      print*,'Error in safely_open:'
      print*,'Iostat=',iost,' while trying to find file ',file_name,'.'
      print*,''
    endif

  endif
!
!
!
! Open file in unit
!--------------------------------------------------------------------!
  if (from_start) then
    open(unit=file_unit,file=file_name,iostat=iost)
  else
    open(unit=file_unit,file=file_name,iostat=iost,access='APPEND')
  endif

  if (iost.ne.0) then
    print*,'Error in safely_open:'
    print*,'Iostat=',iost,' while trying to open file ',file_name,'.'
    if (from_start) then
      print*,'(access attempted: from start)'
    else
      print*,'(access attempted: appending)'
    endif
    print*,''
  endif
!
!
!
  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
