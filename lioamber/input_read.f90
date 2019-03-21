! * read_options     (reads option inputfile.)                                 !
! * read_coords      (reads coordinates inputfile.)                            !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% READ_OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Reads LIO options from an input file.                                        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_options(inputFile, extern_stat)
    use field_subs  , only: read_fields
    use lionml_subs , only: lionml_read, lionml_write
    use garcha_mod  , only: DIIS, hybrid_converg
    use converger_ls, only: Rho_LS

    implicit none
    character(len=20), intent(in)    :: inputFile
    integer, optional, intent(inout) :: extern_stat

    integer :: ios, intern_stat
    logical :: fileExists

    extern_stat = 0
    inquire(file = inputFile, exist = fileExists)
    if(fileExists) then
       open(unit = 100, file = inputFile, iostat = ios)
       call lionml_read(100, intern_stat)
       close(unit = 100)
       extern_stat = intern_stat
       if (intern_stat > 1) return
       call lionml_write
    else
       write(*,*) 'File ', trim(inputFile), ' not found.'
       extern_stat = -3
       return
    endif

    if (Rho_LS .gt. 0 .and. (DIIS .or. hybrid_converg)) then
      hybrid_converg=.false.
      DIIS=.false.
      write(*,*) 'WARNING - read_options: turning off DIIS because of rho_linsearch.'
    end if

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

    integer :: ios, whitespace_count, word_count
    logical :: fileExists
    character(len=1)  :: char_read

    inquire(file=inputCoord,exist=fileExists)
    if(fileExists) then
        open(unit=101,file=inputCoord,iostat=ios)
    else
        write(*,*) 'Input coordinates file ', trim(inputCoord),' not found.'
        stop
    endif

    ! Reads coordinates file.
    ntatom = natom + nsol
    allocate (iz(natom), r(ntatom,3), rqm(natom,3), pc(ntatom))

    ! This is to compatibilize old LIO formats with proper XYZ files.
    whitespace_count = 0
    word_count       = 0
    do while (.true.)
        read(101, '(A1)', advance='no', iostat=ios) char_read
        if (ios /= 0) exit
        if (char_read == " ") then
            whitespace_count = whitespace_count +1
        else if (whitespace_count > 0) then
            word_count       = word_count +1
            whitespace_count = 0
        endif
    enddo
    if ( word_count > 2 ) then
        rewind(101)
    else
        read(101,*)
    endif

    ! Finally reads atomic number / partial charge and coordinates.
    do i=1,natom
        read(101,*) iz(i), r(i,1:3)
        rqm(i,1:3) = r(i,1:3)
    enddo
    do i=natom+1,ntatom
        read(101,*) pc(i), r(i,1:3)
    enddo
    r  = r   / 0.529177D0
    rqm= rqm / 0.529177D0

    close(101)
end subroutine read_coords
