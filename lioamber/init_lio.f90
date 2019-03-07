!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% LIO_INIT.F90  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains common initialization routines, including assignation of  !
! default values for options, wether LIO is run alone or in tantem with AMBER  !
! or GROMACS software packages. Routines currently included are:               !
! * lio_defaults     (called from the last two routines and liomd/liosolo)     !
! * init_lio_common  (called from the last two routines and liomd/liosolo)     !
! * init_lio_gromacs (calls the first two routines when running with GROMACS)  !
! * init_lio_amber   (calls the first two routines when running with AMBER)    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% LIO_DEFAULTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Subroutine lio_defaults gives default values to LIO runtime options.         !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine lio_defaults()

    use garcha_mod, only : fmulliken, fcoord, OPEN, NMAX, DIIS, ndiis, VCINP,  &
                           GOLD, told, Etold, hybrid_converg, good_cut, Iexch, &
                           restart_freq, frestartin, IGRID, frestart, predcoef,&
                           cubegen_only, cube_res, cube_dens, cube_orb,        &
                           cube_sel, cube_orb_file, cube_dens_file, NUNP,      &
                           energy_freq, writeforces, charge, sol, primera,     &
                           cube_elec, cube_elec_file, cube_sqrt_orb, MEMO,     &
                           watermod, fukui, little_cube_size, sphere_radius,   &
                           max_function_exponent, min_points_per_cube,         &
                           assign_all_functions, remove_zero_weights,          &
                           energy_all_iterations, free_global_memory, dipole,  &
                           lowdin, mulliken, print_coeffs, number_restr, Dbug, &
                           steep, Force_cut, Energy_cut, minimzation_steep,    &
                           n_min_steeps, lineal_search, n_points, timers,      &
                           writexyz, IGRID2, propagator, NBCH

    use ECP_mod   , only : ecpmode, ecptypes, tipeECP, ZlistECP, cutECP,       &
                           local_nonlocal, ecp_debug, ecp_full_range_int,      &
                           verbose_ECP, Cnorm, FOCK_ECP_read, FOCK_ECP_write,  &
                           Fulltimer_ECP, cut2_0, cut3_0
    implicit none

!   Names of files used for input and output.
    fmulliken      = 'mulliken'    ; fcoord             = 'qm.xyz'      ;

!   Theory level options.
    OPEN           = .false.       ; told               = 1.0D-6        ;
    NMAX           = 100           ; Etold              = 1.0d0         ;
    good_cut       = 1.0D-3        ; DIIS               = .true.        ;
    ndiis          = 30            ; GOLD               = 10.0D0        ;
    charge         = 0             ; hybrid_converg     = .false.       ;

!   Effective Core Potential options.
    ecpmode        = .false.       ; cut2_0             = 15.d0         ;
    ecptypes       = 0             ; cut3_0             = 12.d0         ;
    tipeECP        = 'NOT-DEFINED' ; verbose_ECP        = 0             ;
    ZlistECP       = 0             ; ecp_debug          = .false.       ;
    FOCK_ECP_read  = .false.       ; Fulltimer_ECP      = .false.       ;
    FOCK_ECP_write = .false.       ; local_nonlocal     = 0             ;
    cutECP         = .true.        ; ecp_full_range_int = .false.       ;

!   TD-DFT options.
    propagator     = 1             ; NBCH               = 10            ;

!   Distance restrain options
    number_restr   = 0             ;

!   Geometry Optimizations
    steep= .false.                 ; Force_cut=1D-5                     ;
    Energy_cut= 1D-4               ; minimzation_steep=5D-2             ;
    n_min_steeps = 500             ; lineal_search=.true.               ;
    n_points = 5                   ;

!   Debug
    Dbug = .false.                 ;

!   Write options and Restart options.
    writexyz       = .true.        ;
    print_coeffs   = .false.       ; frestart           ='restart.out'  ;
    VCINP          = .false.       ; frestartin         = 'restart.in'  ;
    restart_freq   = 0             ; writeforces        = .false.       ;
    fukui          = .false.       ; lowdin             = .false.       ;
    mulliken       = .false.       ; dipole             = .false.       ;

!   Old GPU_options
    max_function_exponent = 10     ; little_cube_size     = 8.0         ;
    min_points_per_cube   = 1      ; assign_all_functions = .false.     ;
    sphere_radius         = 0.6    ; remove_zero_weights  = .true.      ;
    energy_all_iterations = .false.; free_global_memory   = 0.0         ;

!   Cube, grid and other options.
    predcoef       = .false.       ; cubegen_only       = .false.       ;
    cube_res       = 40            ;
    cube_dens      = .false.       ; cube_orb           = .false.       ;
    Iexch          = 9             ; cube_sel           = 0             ;
    cube_orb_file  = "orb.cube"    ; cube_dens_file     = 'dens.cube'   ;
    IGRID          = 2             ; cube_elec          = .false.       ;
    IGRID2         = 2             ; cube_elec_file     = 'field.cube'  ;
    timers         = 0             ; 
    NUNP           = 0             ; energy_freq        = 1             ;
    cube_sqrt_orb  = .false.       ; MEMO               = .true.        ;
    sol            = .false.       ;
    primera        = .true.        ; watermod           = 0             ;

    return
end subroutine lio_defaults
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% INIT_LIO_COMMON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs LIO variable initialization.                                        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine init_lio_common(natomin, Izin, nclatom, callfrom)

    use garcha_mod, only : nunp, RMM, d, r, v, rqm, Em, Rm, pc, Iz, natom, ng0,&
                           ngd0, ngrid, norbit, ntatom, free_global_memory,    &
                           assign_all_functions, energy_all_iterations,        &
                           remove_zero_weights, min_points_per_cube,           &
                           max_function_exponent, little_cube_size,            &
                           OPEN, timers, MO_coef_at, MO_coef_at_b, charge,     &
                           Fmat_vec, Fmat_vec2, Pmat_vec, Hmat_vec, Ginv_vec,  &
                           Gmat_vec, sphere_radius
    use ECP_mod,    only : Cnorm, ecpmode
    use field_data, only : chrg_sq
    use fileio    , only : lio_logo
    use fileio_data, only: style, verbose
    use basis_data, only: M, Md, basis_set, fitting_set, MM, MMd
    use basis_subs, only: basis_init

    implicit none
    integer , intent(in) :: nclatom, natomin, Izin(natomin), callfrom
    integer              :: ng2, ngdDyn, ngDyn, ierr, ios, iostat

    if (verbose .gt. 2) then
      write(*,*)
      write(*,'(A)') "LIO initialisation."
    endif
    call g2g_timer_start('lio_init')
    call g2g_set_options(free_global_memory, little_cube_size, sphere_radius, &
                         assign_all_functions, energy_all_iterations,         &
                         remove_zero_weights, min_points_per_cube,            &
                         max_function_exponent, timers, verbose)

    chrg_sq = charge**2
    if (callfrom.eq.1) then
        natom  = natomin
        if (.not.(allocated(Iz))) allocate(Iz(natom))
        Iz = Izin
        ntatom = natom + nclatom
        allocate(r(ntatom,3), rqm(natom,3), pc(ntatom))
    endif

    ! ngDyn : n° of atoms times the n° of basis functions.                   !
    ! norbit: n° of molecular orbitals involved.                               !
    ! ngdDyn: n° of atoms times the n° of auxiliary functions.               !
    ! Ngrid : n° of grid points (LS-SCF part).                                 !
    ! NOTES: Ngrid may be set to 0  in the case of Numerical Integration. For  !
    ! large systems, ng2 may result in <0 due to overflow.                     !
    call basis_init(basis_set, fitting_set, natom, Iz, iostat)
    if (iostat .gt. 0) then
      stop
      return
   endif

    ngDyn = M; ngdDyn = Md
    ng2 = 5*ngDyn*(ngDyn+1)/2 + 3*ngdDyn*(ngdDyn+1)/2 + &
          ngDyn  + ngDyn*norbit + Ngrid

    allocate(RMM(ng2) , d(natom, natom), v(ntatom,3), Em(ntatom), Rm(ntatom))
    allocate(Fmat_vec(MM), Fmat_vec2(MM), Pmat_vec(MM), Hmat_vec(MM), &
             Ginv_vec(MMd), Gmat_vec(MMd))
    ! Cnorm contains normalized coefficients of basis functions.
    ! Differentiate C for x^2,y^2,z^2 and  xy,xz,yx (3^0.5 factor)
    if (ecpmode) allocate (Cnorm(ngDyn,13))

    call g2g_init()
    allocate(MO_coef_at(ngDyn*ngDyn))
    if (OPEN) allocate(MO_coef_at_b(ngDyn*ngDyn))

    ! Prints chosen options to output.
    call drive(ng2, ngDyn, ngdDyn, iostat)
    if (iostat .gt. 0) then
      stop
      return
    endif

    call g2g_timer_stop('lio_init')

    return
end subroutine init_lio_common
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% INIT_LIO_GROMACS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Subroutine init_lio_gromacs performs Lio initialization when called from     !
! GROMACS software package, in order to conduct a hybrid QM/MM calculation.    !
! In order to avoid compatibility problems due to a FORTRAN/C++ interface, LIO !
! options are read from a file named "lio.in" in the current workspace.        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine init_lio_gromacs(natomin, Izin, nclatom, chargein)
    use garcha_mod, only: charge

    implicit none
    integer,  intent(in) :: chargein, nclatom, natomin, Izin(natomin)
    integer              :: ierr
    character(len=20)    :: inputFile

    ! Gives default values to runtime variables.
    call lio_defaults()
    charge = chargein

    ! Checks if input file exists and writes data to namelist variables.
    inputFile = 'lio.in'
    call read_options(inputFile, ierr)
    if (ierr > 0) return

    ! Initializes LIO. The last argument indicates LIO is not being used alone.
    call init_lio_common(natomin, Izin, nclatom, 1)
end subroutine init_lio_gromacs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% INIT_LIO_AMBER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Subroutine init_lio_amber performs Lio initialization when called from AMBER !
! software package, in order to conduct a hybrid QM/MM calculation.            !
! AMBER directly passes options to LIO, but since the interface has not been   !
! officialy updated on the AMBER side, only some variables are received.       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine init_lio_amber(natomin, Izin, nclatom, charge_i, basis_i            &
           , output_i, fcoord_i, fmulliken_i, frestart_i, frestartin_i         &
           , verbose_i, OPEN_i, NMAX_i, NUNP_i, VCINP_i, GOLD_i, told_i        &
           , rmax_i, rmaxs_i, predcoef_i, idip_i, writexyz_i                   &
           , intsoldouble_i, DIIS_i, ndiis_i, dgtrig_i, Iexch_i, integ_i       &
           , DENS_i , IGRID_i, IGRID2_i , timedep_i , tdstep_i                 &
           , ntdstep_i, field_i, exter_i, a0_i, epsilon_i, Fx_i, Fy_i          &
           , Fz_i, NBCH_i, propagator_i, writedens_i, tdrestart_i              &
           )

    use garcha_mod , only: fmulliken, fcoord, OPEN, NMAX, charge, DIIS,ndiis,&
                           GOLD, told, Etold, hybrid_converg, good_cut,      &
                           propagator, NBCH, VCINP, restart_freq, writexyz,  &
                           frestart, predcoef, frestartin, energy_freq,      &
                           IGRID, IGRID2, nunp, iexch
    use td_data    , only: tdrestart, tdstep, ntdstep, timedep, writedens
    use field_data , only: field, a0, epsilon, Fx, Fy, Fz
    use basis_data , only: int_basis, rmax, rmaxs, basis_set
    use fileio_data, only: verbose
    use ECP_mod    , only: ecpmode, ecptypes, tipeECP, ZlistECP, cutECP,     &
                           local_nonlocal, ecp_debug, ecp_full_range_int,    &
                           verbose_ECP, Cnorm, FOCK_ECP_read, FOCK_ECP_write,&
                           Fulltimer_ECP, cut2_0, cut3_0

    implicit none
    integer , intent(in) :: charge_i, nclatom, natomin, Izin(natomin)
    character(len=20) :: basis_i, fcoord_i, fmulliken_i, frestart_i, &
                         frestartin_i, inputFile
    logical           :: verbose_i, OPEN_i, VCINP_i, predcoef_i, writexyz_i,   &
                         DIIS_i, field_i, exter_i, writedens_i, tdrestart_i
    integer           :: NMAX_i, NUNP_i, ndiis_i, Iexch_i, IGRID_i, IGRID2_i,  &
                         timedep_i, ntdstep_i, NBCH_i, propagator_i, ierr
    real*8            :: GOLD_i, told_i, rmax_i, rmaxs_i, tdstep_i,  &
                         a0_i, epsilon_i, Fx_i, Fy_i, Fz_i
    ! Deprecated or removed variables
    character(len=20) :: output_i
    integer           :: idip_i
    logical           :: intsoldouble_i, dens_i, integ_i
    double precision  :: dgtrig_i

    ! Gives default values to variables.
    call lio_defaults()

    fcoord         = fcoord_i       ; fmulliken     = fmulliken_i    ;
    frestart       = frestart_i     ; frestartin    = frestartin_i   ;
    OPEN           = OPEN_i         ;
    NMAX           = NMAX_i         ; NUNP          = NUNP_i         ;
    VCINP          = VCINP_i        ; GOLD          = GOLD_i         ;
    told           = told_i         ; rmax          = rmax_i         ;
    rmaxs          = rmaxs_i        ; predcoef      = predcoef_i     ;
    writexyz       = writexyz_i     ;
    DIIS           = DIIS_i         ; ndiis         = ndiis_i        ;
    Iexch          = Iexch_i        ;
    IGRID          = IGRID_i        ;
    IGRID2         = IGRID2_i       ; timedep       = timedep_i      ;
    field          = field_i        ; tdrestart     = tdrestart_i    ;
    tdstep         = tdstep_i       ; ntdstep       = ntdstep_i      ;
    a0             = a0_i           ; epsilon       = epsilon_i      ;
    Fx             = Fx_i           ; Fy            = Fy_i           ;
    Fz             = Fz_i           ; NBCH          = NBCH_i         ;
    propagator     = propagator_i   ; writedens     = writedens_i    ;
    charge         = charge_i       ;
    if (verbose_i) verbose = 1

    ! Checks if input file exists and writes data to namelist variables. 
    ! Previous options are overwritten.
    inputFile = 'lio.in'
    call read_options(inputFile, ierr)
    if (ierr > 0) return
    if ((.not. int_basis) .and. (basis_i .ne. 'basis')) basis_set = basis_i


    ! Initializes LIO. The last argument indicates LIO is not being used alone.
    call init_lio_common(natomin, Izin, nclatom, 1)

    return
end subroutine init_lio_amber
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% INIT_LIOAMBER_EHREN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs LIO variable initialization.                                        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine init_lioamber_ehren(natomin, Izin, nclatom, charge_i, basis_i       &
           , output_i, fcoord_i, fmulliken_i, frestart_i, frestartin_i         &
           , verbose_i, OPEN_i, NMAX_i, NUNP_i, VCINP_i, GOLD_i, told_i        &
           , rmax_i, rmaxs_i, predcoef_i, idip_i, writexyz_i                   &
           , intsoldouble_i, DIIS_i, ndiis_i, dgtrig_i, Iexch_i, integ_i       &
           , DENS_i , IGRID_i, IGRID2_i , timedep_i , tdstep_i                 &
           , ntdstep_i, field_i, exter_i, a0_i, epsilon_i, Fx_i, Fy_i          &
           , Fz_i, NBCH_i, propagator_i, writedens_i, tdrestart_i, dt_i        &
           )

   use garcha_mod , only: first_step, doing_ehrenfest
   use basis_subs , only: basis_setup_ehren
   use td_data    , only: timedep, tdstep
   use lionml_data, only: ndyn_steps, edyn_steps
   use liosubs    , only: catch_error


   implicit none
   integer, intent(in) :: charge_i, nclatom, natomin, Izin(natomin)

   character(len=20) :: basis_i, output_i, fcoord_i, fmulliken_i, frestart_i   &
                     &, frestartin_i, inputFile

   logical :: verbose_i, OPEN_i, VCINP_i, predcoef_i, writexyz_i, DIIS_i       &
           &, intsoldouble_i, integ_i, DENS_i, field_i, exter_i, writedens_i   &
           &, tdrestart_i

   integer :: NMAX_i, NUNP_i, idip_i, ndiis_i, Iexch_i, IGRID_i, IGRID2_i      &
           &, timedep_i, ntdstep_i, NBCH_i, propagator_i

   real*8  :: GOLD_i, told_i, rmax_i, rmaxs_i, dgtrig_i, tdstep_i, a0_i        &
           &, epsilon_i, Fx_i, Fy_i, Fz_i, dt_i


   call init_lio_amber(natomin, Izin, nclatom, charge_i, basis_i               &
           , output_i, fcoord_i, fmulliken_i, frestart_i, frestartin_i         &
           , verbose_i, OPEN_i, NMAX_i, NUNP_i, VCINP_i, GOLD_i, told_i        &
           , rmax_i, rmaxs_i, predcoef_i, idip_i, writexyz_i                   &
           , intsoldouble_i, DIIS_i, ndiis_i, dgtrig_i, Iexch_i, integ_i       &
           , DENS_i , IGRID_i, IGRID2_i , timedep_i , tdstep_i                 &
           , ntdstep_i, field_i, exter_i, a0_i, epsilon_i, Fx_i, Fy_i          &
           , Fz_i, NBCH_i, propagator_i, writedens_i, tdrestart_i              &
           )
   call basis_setup_ehren()

   first_step=.true.

   if ( (ndyn_steps>0) .and. (edyn_steps>0) ) doing_ehrenfest=.true.

   tdstep = (dt_i) * (41341.3733366d0)

!  Amber should have time units in 1/20.455 ps, but apparently it has time
!  in ps. Just have to transform to atomic units
!  ( AU = 2.418884326505 x 10e-17 s )

end subroutine init_lioamber_ehren
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% INIT_LIO_HYBRID  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Subroutine init_lio_hybrid performs Lio initialization when called from      !
! Hybrid software package, in order to conduct a hybrid QM/MM calculation.     !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine init_lio_hybrid(hyb_natom, mm_natom, chargein, iza)
    use garcha_mod, only: charge

    implicit none
    integer, intent(in) :: hyb_natom !number of total atoms
    integer, intent(in) :: mm_natom  !number of MM atoms
    integer             :: ierr
    character(len=20)   :: inputFile
    integer, intent(in) :: chargein   !total charge of QM system
    integer, dimension(hyb_natom), intent(in) :: iza  !array of charges of all QM/MM atoms

    ! Gives default values to runtime variables.
    call lio_defaults()
    charge = chargein

    ! Checks if input file exists and writes data to namelist variables.
    inputFile = 'lio.in'
    call read_options(inputFile, ierr)
    if (ierr > 0) return

    ! Initializes LIO. The last argument indicates LIO is not being used alone.
    call init_lio_common(hyb_natom, Iza, mm_natom, 1)

    return
end subroutine init_lio_hybrid
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
