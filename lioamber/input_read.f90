! * read_options     (reads option inputfile.)                                 !
! * read_coords      (reads coordinates inputfile.)                            !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% READ_OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Reads LIO options from an input file.                                        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_options(inputFile)
    use field_subs , only: read_fields
    use lionml_subs, only: lionml_read, lionml_write

    implicit none
    character(len=20), intent(in)  :: inputFile

    integer :: ios
    logical :: fileExists

    inquire(file = inputFile, exist = fileExists)
    if(fileExists) then
       open(unit = 100, file = inputFile, iostat = ios)
       call lionml_read( 100 )
       close(unit = 100)
       call lionml_write
    else
       write(*,*) 'File ', adjustl(inputFile), ' not found. Using defaults.'
    endif

    call read_fields()

    return
end subroutine read_options
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% READ_COORDS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Reads atoms' coordinates from an input file.                                 !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_coords(inputCoord)

    use garcha_mod, only : natom, ntatom, nsol, iz, r, rqm, pc

    character(len=20), intent(in) :: inputCoord

    integer :: ios
    logical :: fileExists

    inquire(file=inputCoord,exist=fileExists)
    if(fileExists) then
        open(unit=101,file=inputCoord,iostat=ios)
    else
        write(*,*) 'Input coordinates file ',adjustl(inputCoord),' not found.'
        stop
    endif

    ! Reads coordinates file.
    ntatom = natom + nsol
    allocate (iz(natom), r(ntatom,3), rqm(natom,3), pc(ntatom))
    do i=1,natom
        read(101,*) iz(i), r(i,1:3)
        rqm(i,1:3) = r(i,1:3)
    enddo
    do i=natom+1,ntatom
        read(101,*) pc(i), r(i,1:3)
    enddo
    r  = r   / 0.529177D0
    rqm= rqm / 0.529177D0

    return
end subroutine read_coords
