! * read_options     (reads option inputfile.)                                 !
! * read_coords      (reads coordinates inputfile.)                            !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% READ_OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Reads LIO options from an input file.                                        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_options(inputFile, extern_stat)
    use converger_subs, only: converger_options_check
    use cdft_subs     , only: cdft_options_check, cdft_input_read
    use field_subs    , only: read_fields
    use garcha_mod    , only: energy_all_iterations, becke, open
    use lionml_subs   , only: lionml_read, lionml_write

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
       call cdft_input_read(100)
       close(unit = 100)
       extern_stat = intern_stat
       if (intern_stat > 1) return
    else
       write(*,*) 'File ', trim(inputFile), ' not found.'
       extern_stat = -3
       return
    endif

    call converger_options_check(energy_all_iterations)
    call cdft_options_check(becke, open)
    call lionml_write()
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
    use constants_mod, only : bohr
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
    r  = r   / bohr
    rqm= rqm / bohr
    call recenter_coords(rqm, r, natom, nsol)

    close(101)
end subroutine read_coords

! Takes the geometric center of the QM system and rescales everything.
subroutine recenter_coords(pos_qm, pos_tot, n_qm, n_sol)
    implicit none
    integer     , intent(in)    :: n_qm, n_sol
    real(kind=8), intent(inout) :: pos_qm(n_qm,3), pos_tot(n_qm+n_sol,3)

    real(kind=8) :: geom_center(3) = 0.0D0
    integer      :: iatom

    ! Gets geometric center of the QM system. This way we avoid
    ! recalculating everything if there are MM atoms.
    do iatom = 1, n_qm
        geom_center(1) = geom_center(1) + pos_qm(iatom,1)
        geom_center(2) = geom_center(2) + pos_qm(iatom,2)
        geom_center(3) = geom_center(3) + pos_qm(iatom,3)
    enddo
    geom_center = geom_center / dble(n_qm)
    do iatom = 1, n_qm
        pos_qm(iatom,1)  = pos_qm(iatom,1)  - geom_center(1)
        pos_qm(iatom,2)  = pos_qm(iatom,2)  - geom_center(2)
        pos_qm(iatom,3)  = pos_qm(iatom,3)  - geom_center(3)
        pos_tot(iatom,1) = pos_tot(iatom,1) - geom_center(1)
        pos_tot(iatom,2) = pos_tot(iatom,2) - geom_center(2)
        pos_tot(iatom,3) = pos_tot(iatom,3) - geom_center(3)
    enddo

    if (n_sol < 1) return
    do iatom = n_qm +1, n_qm + n_sol
        pos_tot(iatom,1) = pos_tot(iatom,1) - geom_center(1)
        pos_tot(iatom,2) = pos_tot(iatom,2) - geom_center(2)
        pos_tot(iatom,3) = pos_tot(iatom,3) - geom_center(3)
    enddo 

end subroutine recenter_coords
