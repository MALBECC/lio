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

    use garcha_mod, only : fmulliken, fcoord, OPEN, VCINP,                     &
                           Iexch, restart_freq, frestartin, IGRID, frestart,   &
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
                           calc_propM, writexyz, IGRID2, propagator, NBCH,     &
                           predcoef

    use ECP_mod   , only : ecpmode, ecptypes, tipeECP, ZlistECP, cutECP,       &
                           local_nonlocal, ecp_debug, ecp_full_range_int,      &
                           verbose_ECP, FOCK_ECP_read, FOCK_ECP_write,         &
                           Fulltimer_ECP, cut2_0, cut3_0, first_steep

    use rhoint   , only : write_int_rho, w_rho_xmin, w_rho_ymin, w_rho_zmin,   &
                          w_rho_xmax, w_rho_ymax, w_rho_zmax, w_rho_dx,        &
                          w_rho_dy, w_rho_dz, w_rho_rmin, w_rho_rmax, w_rho_dr,&
                          w_rho_dtheta, w_rho_dphi, write1Drho


    implicit none

!   Names of files used for input and output.
    fmulliken      = 'mulliken'    ; fcoord             = 'qm.xyz'      ;

!   Theory level options.
    OPEN           = .false.       ; charge             = 0             ;

!   Effective Core Potential options.
    ecpmode        = .false.       ; cut2_0             = 35.d0         ;
    ecptypes       = 0             ; cut3_0             = 45.d0         ;
    tipeECP        = 'NOT-DEFINED' ; verbose_ECP        = 0             ;
    ZlistECP       = 0             ; ecp_debug          = .false.       ;
    FOCK_ECP_read  = .false.       ; Fulltimer_ECP      = .false.       ;
    FOCK_ECP_write = .false.       ; local_nonlocal     = 0             ;
    cutECP         = .true.        ; ecp_full_range_int = .false.       ;
    first_steep    = .true.        ;

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
    print_coeffs   = .false.       ; frestart           = 'restart.out' ;
    VCINP          = .false.       ; frestartin         = 'restart.in'  ;
    restart_freq   = 0             ; writeforces        = .false.       ;
    fukui          = .false.       ; lowdin             = .false.       ;
    mulliken       = .false.       ; dipole             = .false.       ;
    print_coeffs   = .false.       ; calc_propM         = .false.       ;
    write_int_rho  = '  '          ; w_rho_xmin         = -5.d0         ;
    w_rho_ymin     = -5.d0         ; w_rho_zmin         = -5.d0         ;
    w_rho_xmax     = 5.d0          ; w_rho_ymax         =  5.d0         ;
    w_rho_zmax     = 5.d0          ; w_rho_dx           =  0.1d0        ;
    w_rho_dy       = 0.1d0         ; w_rho_dz           =  0.1d0        ;
    w_rho_rmin     = 0.d0          ; w_rho_rmax         =  5.d0         ;
    w_rho_dr       = 0.1d0         ; w_rho_dtheta       =  0.1d0        ;
    w_rho_dphi     = 0.1d0         ; write1Drho         = .false.       ;

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

    use garcha_mod, only : d, r, v, rqm, Em, Rm, pc, Iz, natom,                &
                           ntatom, free_global_memory,                         &
                           assign_all_functions, energy_all_iterations,        &
                           remove_zero_weights, min_points_per_cube,           &
                           max_function_exponent, little_cube_size,            &
                           OPEN, timers, MO_coef_at, MO_coef_at_b, charge,     &
                           Fmat_vec, Fmat_vec2, Pmat_vec, Hmat_vec, Ginv_vec,  &
                           Gmat_vec, Pmat_en_wgt, sphere_radius
    use ECP_mod,    only : Cnorm, ecpmode
    use field_data, only : chrg_sq
    use fileio    , only : lio_logo
    use fileio_data, only: verbose
    use basis_data, only: M, basis_set, fitting_set, MM, MMd
    use basis_subs, only: basis_init
    use tbdft_data, only: MTB, tbdft_calc
    use dftd3     , only: dftd3_setup

    implicit none
    integer , intent(in) :: nclatom, natomin, Izin(natomin), callfrom
    integer              :: iostat
    integer              :: M_f

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
    
    ! TBDFT: Updating M for TBDFT calculations
    M_f = M
    if (tbdft_calc /= 0) M_f = M + MTB

    if (iostat > 0) then
      stop
      return
   endif

    allocate(d(natom, natom), v(ntatom,3), Em(ntatom), Rm(ntatom))
    allocate(Fmat_vec(MM), Fmat_vec2(MM), Pmat_vec(MM), Hmat_vec(MM), &
             Ginv_vec(MMd), Gmat_vec(MMd), Pmat_en_wgt(MM))

    d  = 0.0D0 ; v  = 0.0D0 ; Em = 0.0D0 ; Rm = 0.0D0
    Pmat_vec    = 0.0D0 ; Hmat_vec  = 0.0D0 ; Fmat_vec  = 0.0D0
    Fmat_vec2   = 0.0D0 ; Gmat_vec  = 0.0D0 ; Ginv_vec  = 0.0D0
    Pmat_en_wgt = 0.0D0
    ! Cnorm contains normalized coefficients of basis functions.
    ! Differentiate C for x^2,y^2,z^2 and  xy,xz,yx (3^0.5 factor)
    if (ecpmode) allocate (Cnorm(M,13))

    call g2g_init()
    allocate(MO_coef_at(M_f,M_f))
    if (OPEN) allocate(MO_coef_at_b(M_f,M_f))

    ! Prints chosen options to output.
    call drive(iostat)
    if (iostat .gt. 0) then
      stop
      return
    endif
    call dftd3_setup(natom, Iz)

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

    use garcha_mod , only: fmulliken, fcoord, OPEN, charge, propagator, NBCH, &
                           VCINP, writexyz, frestart, predcoef, frestartin,   &
                           IGRID, IGRID2, nunp, iexch
    use td_data    , only: tdrestart, tdstep, ntdstep, timedep, writedens
    use field_data , only: field, a0, epsilon, Fx, Fy, Fz
    use basis_data , only: int_basis, rmax, rmaxs, basis_set
    use fileio_data, only: verbose
    use converger_data, only: DIIS, nDIIS, gOld, tolD, nMax

    implicit none
    ! Variables received from &lio namelist in amber input file.
    character(len=20), intent(in) :: basis_i, fcoord_i, fmulliken_i, &
                                     frestart_i, frestartin_i
    integer          , intent(in) :: charge_i, nclatom, natomin, Izin(natomin),&
                                     NMAX_i, NUNP_i, ndiis_i, Iexch_i, IGRID_i,&
                                     IGRID2_i, timedep_i, ntdstep_i, NBCH_i,   &
                                     propagator_i
    logical          , intent(in) :: verbose_i, OPEN_i, VCINP_i, predcoef_i,   &
                                     writexyz_i, DIIS_i, field_i, exter_i,     &
                                     writedens_i, tdrestart_i
    real(kind=8)     , intent(in) :: GOLD_i, told_i, rmax_i, rmaxs_i, tdstep_i,&
                                     a0_i, epsilon_i, Fx_i, Fy_i, Fz_i

    ! Deprecated or removed variables
    character(len=20), intent(in) :: output_i
    integer          , intent(in) :: idip_i
    logical          , intent(in) :: intsoldouble_i, dens_i, integ_i
    real(kind=8)     , intent(in) :: dgtrig_i

    character(len=20) :: inputFile
    integer           :: ierr
    logical           :: file_exists

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

    inquire(file = basis_i, exist = file_exists)
    if (file_exists) then
        write(*,'(A)') "LIO - Custom basis set found, using present file."
        int_basis = .false.
        basis_set = basis_i
    endif

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
   use td_data    , only: tdstep
   use lionml_data, only: ndyn_steps, edyn_steps
   use liosubs    , only: catch_error
   use excited_data,only: TSH, tsh_time_dt, tsh_coef, tsh_Jstate, &
                          tsh_Kstate, gamma_old, excited_forces

   implicit none
   integer, intent(in) :: charge_i, nclatom, natomin, Izin(natomin)

   character(len=20) :: basis_i, output_i, fcoord_i, fmulliken_i, frestart_i   &
                     &, frestartin_i

   logical :: verbose_i, OPEN_i, VCINP_i, predcoef_i, writexyz_i, DIIS_i       &
           &, intsoldouble_i, integ_i, DENS_i, field_i, exter_i, writedens_i   &
           &, tdrestart_i

   integer :: NMAX_i, NUNP_i, idip_i, ndiis_i, Iexch_i, IGRID_i, IGRID2_i      &
           &, timedep_i, ntdstep_i, NBCH_i, propagator_i

   real*8  :: GOLD_i, told_i, rmax_i, rmaxs_i, dgtrig_i, tdstep_i, a0_i        &
           &, epsilon_i, Fx_i, Fy_i, Fz_i, dt_i

   ! SEED VARIABLES
   integer :: random_size
   integer, dimension(12) :: random_values
   integer, dimension(:), allocatable :: seed

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

   if ( TSH ) then
      ! dt_i = ps
      ! 1 ps = 4.134137d4 au
      ! tsh_time_dt = au

      ! RANDOM SEED
      call date_and_time(VALUES=random_values)
      call random_seed(size=random_size)
      allocate(seed(random_size))
      seed = random_values
      print*, "SEED:",seed
      call random_seed(put=seed)
      deallocate(seed)

      print*, "*Init TSH Dynamics"
      tsh_time_dt = tdstep ! tsh_time_dt in atomic units
      allocate(tsh_coef(2))
      tsh_coef(1) = (0.0d0,0.0d0)
      tsh_coef(2) = (1.0d0,0.0d0)
      tsh_Jstate  = 2
      tsh_Kstate  = 1
      allocate(gamma_old(natomin,3))
      gamma_old = 0.0d0

      excited_forces = .true.
   endif

end subroutine init_lioamber_ehren
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% INIT_LIO_HYBRID  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Subroutine init_lio_hybrid performs Lio initialization when called from      !
! Hybrid software package, in order to conduct a hybrid QM/MM calculation.     !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine init_lio_hybrid(version_check, hyb_natom, mm_natom, chargein, iza, spin)
    use garcha_mod, only: OPEN, Nunp, charge

    implicit none
    integer, intent(in) :: hyb_natom !number of total atoms
    integer, intent(in) :: mm_natom  !number of MM atoms
    integer, intent(in) :: version_check !check version compatibility
    integer             :: ierr
    character(len=20)   :: inputFile
    integer, intent(in) :: chargein   !total charge of QM system
    integer, dimension(hyb_natom), intent(in) :: iza  !array of charges of all QM/MM atoms
    double precision, intent(in) :: spin !number of unpaired electrons
    integer :: Nunp_aux !auxiliar

    if (version_check.ne.1) Stop 'LIO version is not compatible with hybrid'

    ! Gives default values to runtime variables.
    call lio_defaults()
    charge = chargein

    !select spin case
    Nunp_aux=int(spin)
    Nunp=Nunp_aux

    ! Checks if input file exists and writes data to namelist variables.
    inputFile = 'lio.in'
    call read_options(inputFile, ierr)
    if (ierr > 0) return
    !select spin case
    Nunp_aux=int(spin)
    if (Nunp_aux .ne. Nunp) STOP "lio.in have a different spin than *.fdf"
    if (Nunp .ne. 0) OPEN=.true.
    if (OPEN) write(*,*) "Runing hybrid open shell, with ", Nunp, "unpaired electrons"

    if (Nunp_aux .ne. Nunp) STOP "lio.in have a different spin than *.fdf"
    if (Nunp .ne. 0) OPEN=.true.
    if (OPEN) write(*,*) "Runing hybrid open shell, with ", Nunp, "unpaired electrons"

    ! Initializes LIO. The last argument indicates LIO is not being used alone.
    call init_lio_common(hyb_natom, Iza, mm_natom, 1)

    return
end subroutine init_lio_hybrid
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
