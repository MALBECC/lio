! * read_options     (reads option inputfile.)                                 !
! * read_coords      (reads coordinates inputfile.)                            !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% READ_OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Reads LIO options from an input file.                                        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_options(inputFile, charge)

    use garcha_mod, only : natom, nsol, basis, output, fmulliken, fcoord, OPEN,&
                           NMAX, basis_set, fitting_set, int_basis, DIIS,      &
                           ndiis, GOLD, told, Etold, hybrid_converg, good_cut, &
                           rmax, rmaxs, omit_bas, propagator, NBCH, Fx, Fy, Fz,&
                           field, verbose, VCINP, restart_freq, writexyz,      &
                           frestartin, frestart, predcoef, idip, intsoldouble, &
                           dgtrig, Iexch, integ, DENS, IGRID, IGRID2, epsilon, &
                           a0, cubegen_only, cube_res, cube_dens, cube_orb,    &
                           cube_sel, cube_orb_file, cube_dens_file, cube_elec, &
                           cube_elec_file, energy_freq, NUNP, style, allnml,   &
                           writeforces, cube_sqrt_orb, fukui, little_cube_size,&
                           max_function_exponent, min_points_per_cube,         &
                           assign_all_functions, remove_zero_weights,          &
                           energy_all_iterations, free_global_memory,          &
                           sphere_radius, dipole, lowdin, mulliken,            &
                           print_coeffs, number_restr, Dbug, steep, Force_cut, &
                           Energy_cut, minimzation_steep, n_min_steeps,        &
                           lineal_search, n_points, timers
   use td_data    , only : tdrestart, writedens, td_rst_freq, tdstep, ntdstep, &
                           timedep
    use ECP_mod   , only : ecpmode, ecptypes, tipeECP, ZlistECP, verbose_ECP,  &
                           cutECP, local_nonlocal, ecp_debug, FOCK_ECP_read,   &
                           FOCK_ECP_write, ecp_full_range_int, Fulltimer_ECP,  &
                           cut2_0, cut3_0


    use transport_data, only  : transport_calc, generate_rho0, gate_field,     &
                                save_charge_freq, driving_rate, Pop_Drive
    use dftb_data ,only: dftb_calc, MTB, alfaTB, betaTB, gammaTB, Vbias,      &
                         end_basis, start_tdtb, end_tdtb, TBsave, TBload
    implicit none
    character(len=20), intent(in)  :: inputFile
    integer          , intent(out) :: charge

    integer :: ios, iErr
    logical :: fileExists

                   ! Common LIO variables.
    namelist /lio/ OPEN, NMAX, Nunp, VCINP, GOLD, told, Etold, rmax, rmaxs,    &
                   predcoef, idip, writexyz, intsoldouble, DIIS, ndiis, dgtrig,&
                   Iexch, integ, dens, igrid, igrid2, good_cut, hybrid_converg,&
                   ! File Input/Output.
                   frestartin, style, allnml, frestart, fukui, dipole, lowdin, &
                   mulliken, writeforces, int_basis, fitting_set, basis_set,   &
                   restart_freq, print_coeffs,                                 &
                   ! DFT and TD-DFT Variables.
                   timedep, tdstep, ntdstep, propagator, NBCH, field, epsilon, &
                   a0, Fx, Fy, Fz, tdrestart, writedens, td_rst_freq,          &
                   ! Effective Core Potential Variables.
                   ecpmode, ecptypes, tipeECP, ZlistECP, cutECP, ecp_debug,    &
                   local_nonlocal, ecp_debug, ecp_full_range_int, verbose_ECP, &
                   verbose, FOCK_ECP_read, FOCK_ECP_write, Fulltimer_ECP,      &
                   cut2_0, cut3_0,                                             &
                   ! Distance Restrain
                   number_restr,                                               &
                   ! Debug variables
                   Dbug, timers,                                               &
                   ! Geometry optimizations
                   steep, Force_cut, Energy_cut, minimzation_steep,            &
                   n_min_steeps,lineal_search,n_points,                        &
                   ! Variables for orbital printing.
                   cubegen_only, cube_res, cube_sel, cube_dens, cube_dens_file,&
                   cube_orb, cube_orb_file, cube_elec, cube_elec_file,         &
                   cube_sqrt_orb,                                              &
                   ! Variables for GPU options.
                   little_cube_size, max_function_exponent, free_global_memory,&
                   min_points_per_cube, assign_all_functions, sphere_radius,   &
                   remove_zero_weights, energy_all_iterations,                 &
                   ! Variables when LIO is used alone.
                   natom, nsol, charge,                                        &
                   ! Variables for Transport
                   transport_calc, generate_rho0, gate_field,                  &
                   save_charge_freq, driving_rate, Pop_Drive,                  &
                   !Variables for DFTB
                   dftb_calc, MTB, alfaTB, betaTB, gammaTB, Vbias, end_basis,  &
                   start_tdtb, end_tdtb, TBsave, TBload

    inquire(file = inputFile, exist = fileExists)
    if(fileExists) then
        open(unit = 100, file = inputFile, iostat = ios)
        read(100, nml = lio, iostat = iErr)
        if(ierr.gt.0) stop 'Input error in LIO namelist.'
        close(unit = 100)
    else
        write(*,*) 'File ', adjustl(inputFile), ' not found. Using defaults.'
    endif

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
