!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
program liomd
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Make an NVE ensamble MD.

    use garcha_mod, only : basis, basis_set, fitting_set, natom, ntatom, nsol, &
                           Iz, r, rqm, int_basis, omit_bas
    use liosubs   , only : write_energy, write_geom, set_masses, nuclear_verlet
    use fileio_data, only: verbose
    use ehrensubs
    use basis_data

    implicit none
    real*8  :: escf, dipxyz(3), timeStep, Kenergy, Uenergy
    integer :: ii, kk, narg, nqnuc, ios, nn, md_steps, charge
    logical :: filexist

    character(len=20)   :: argument, inpfile, inpbasis, inpcoords
    real*8, allocatable :: oldpos(:,:), nowpos(:,:), newpos(:,:),  &
                           nowvel(:,:), nowfrc(:,:), atom_mass(:), &
                           dxyzqm(:,:), dxyzcl(:,:)

    ! Loads default variable values.
    call lio_defaults()

    ! Reads command line arguments for LIO.
    narg = command_argument_count()
    do ii = 1, narg
        call get_command_argument(ii, argument)
        select case(adjustl(argument))
            case("-i")
                call get_command_argument(ii + 1, inpfile)
            case("-b")
                call get_command_argument(ii + 1, basis)
            case("-bs")
                omit_bas = .true.
                call get_command_argument(ii + 1, basis_set)
            case("-fs")
                call get_command_argument(ii + 1, fitting_set)
            case("-ib")
                int_basis = .true.
            case("-c")
                call get_command_argument(ii + 1, inpcoords)
            case("-v")
                verbose = 4
           case default
        endselect
    enddo

    ! Reads runtime options and input coordinates.
    call read_options(inpfile, charge)
    call read_coords(inpcoords)

    ! Initializes variables.
    call init_lio_common(natom, Iz, nsol, charge, 0)

    ! ALLOCATION AND INITIALIZATION
    if (allocated(dxyzqm)) deallocate (dxyzqm)
    allocate(dxyzqm(3, natom), oldpos(natom, 3), nowpos(natom, 3), &
             newpos(natom, 3), nowvel(natom, 3), nowfrc(natom, 3))

    if (allocated(atom_mass)) deallocate(atom_mass)
    allocate(atom_mass(natom))
    call set_masses(natom, Iz, atom_mass)

    ! Setup. Timestep should be in ps.
    do kk = 1, 3
        do nn = 1, natom
            oldpos(nn, kk)= r(nn, kk)
            nowpos(nn, kk)= r(nn, kk)
            newpos(nn, kk)= r(nn, kk)
            nowvel(nn, kk)= 0.0d0
            nowfrc(nn, kk)= 0.0d0
        enddo
    enddo

    md_steps  = 0
    timeStep = 0.0001                 ! In ps.
    timeStep = timeStep * 41341.0d0  ! In A.U.
    open(unit=501, file='liomd-trays.xyz' )
    open(unit=502, file='liomd-energy.dat')
    call write_energy(-1.0d0, escf, Kenergy, escf + Kenergy, 502)

    do nn = 1, md_steps + 1
        call liomain(escf, dipxyz)
        Uenergy = escf
        call dft_get_qm_forces(dxyzqm)
        do kk = 1, 3
            do ii= 1, natom
                nowfrc(ii, kk)= -dxyzqm(kk, ii)
            enddo
        enddo
        call nuclear_verlet(natom, timeStep, atom_mass, nowfrc,   &
                            oldpos, nowpos, newpos, nowvel, Kenergy)
        call write_geom(natom, Iz, nowpos, 501)
        call write_energy(nn*timeStep, escf, Kenergy, escf+Kenergy, 502)

        oldpos = nowpos
        nowpos = newpos
        do kk = 1, 3
            do ii = 1, natom
                r(ii, kk)   = nowpos(ii, kk)
                rqm(ii, kk) = nowpos(ii, kk)
            enddo
        enddo
    enddo

    write(*,*) 'SCF ENRGY=', escf

    call lio_finalize()

100 format (I5, 2x, f10.6, 2x, f10.6, 2x, f10.6)
end program
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
