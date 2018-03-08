! MAIN SUBROUTINE ----------------------------------------------------
! DFT calculation with gaussian basis sets
!---------------------------------------------------------------------
program liosolo

    use garcha_mod , only : natom, nsol, Iz, basis, basis_set, fitting_set, &
                             int_basis, omit_bas, verbose, writeforces
    use ECP_mod    , only : ecpmode
#ifdef CUBLAS
    use cublasmath
#endif

    implicit none 
    character(len=20) :: argument, inpfile, inpbasis, inpcoords
    integer           :: charge, ifind, ierr, i, k, narg, ios
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
                call get_command_argument(i+1,basis)
            case("-bs")
                omit_bas=.true.
                call get_command_argument(i+1,basis_set)
            case("-fs")
                call get_command_argument(i+1,fitting_set)
            case("-ib")
                int_basis=.true.
            case("-c")
                call get_command_argument(i+1,inpcoords)
            case("-v")
                verbose=.true.
            case default
        endselect
    enddo

    ! Reads options and coordinates files.
    call read_options(inpfile, charge)

    call read_coords(inpcoords)
     
    ! Initializes LIO. The last argument indicates LIO is being used alone.
    call init_lio_common(natom, Iz, nsol, charge, 0)

    ! Calls main procedures.
    call liomain(escf, dipxyz)
    call lio_finalize()

end program liosolo

