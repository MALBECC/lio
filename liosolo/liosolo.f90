! MAIN SUBROUTINE ----------------------------------------------------
! DFT calculation with gaussian basis sets
!---------------------------------------------------------------------
program liosolo

    use garcha_mod , only : natom, nsol, Iz, writeforces, writexyz
    use basis_data , only : basis_set, fitting_set, int_basis
    use ECP_mod    , only : ecpmode
    use fileio_data, only : verbose
    use fileio     , only : lio_logo

    implicit none
    character(len=20) :: argument, inpfile, inpbasis, inpcoords
    integer           :: ifind, ierr, i, k, narg, ios
    real*8            :: dipxyz(3), escf

    ! Calls default values for variables.
    call lio_defaults()

    ! Reads command line arguments for LIO.
    narg=command_argument_count()
    do i=1, narg
        call get_command_argument(i,argument)
        select case(adjustl(argument))
            case("-i")
                call get_command_argument(i+1,inpfile)
            case("-b")
                call get_command_argument(i+1,basis_set)
            case("-bs")
                call get_command_argument(i+1,basis_set)
            case("-fs")
                call get_command_argument(i+1,fitting_set)
            case("-ib")
                int_basis=.true.
            case("-c")
                call get_command_argument(i+1,inpcoords)
            case("-v")
                verbose = 4
            case default
        endselect
    enddo

    ! Prints LIO welcome message.
    call lio_logo()

    ! Reads options and coordinates files. Previously sets writexyz to false in
    ! order to avoid the qm.xyz file when not running MD.
    writexyz = .false.
    ierr = 0
    call read_options(inpfile, ierr)
    select case (ierr)
        case (-1,-4)
            write(*,'(A)')
            write(*,'(A)') "ERROR: &lio namelist is required for LIO-only runs."
            return
        case (-3)
            write(*,'(A)') "ERROR: input file not found."
            return
        case (1:) ! If .ge.1, there was an error in namelists.
            return
        case default
    end select

    call read_coords(inpcoords)

    ! Initializes LIO. The last argument indicates LIO is being used alone.
    call init_lio_common(natom, Iz, nsol, 0)

    ! Calls main procedures.
    call liomain(escf, dipxyz)
    call lio_finalize()

end program liosolo
