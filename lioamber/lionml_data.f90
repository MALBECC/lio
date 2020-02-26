!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
module lionml_data

   use garcha_mod        , only: natom, nsol, fmulliken, fcoord, OPEN,         &
                                 propagator, VCINP, restart_freq, writexyz,    &
                                 Iexch, frestartin, frestart, predcoef,        &
                                 cubegen_only, cube_res, cube_dens, cube_orb,  &
                                 cube_sel, cube_orb_file, cube_dens_file,      &
                                 cube_elec, cube_elec_file, energy_freq, NUNP, &
                                 writeforces, cube_sqrt_orb, NBCH,             &
                                 fukui, little_cube_size, min_points_per_cube, &
                                 max_function_exponent, assign_all_functions,  &
                                 remove_zero_weights, energy_all_iterations,   &
                                 free_global_memory, sphere_radius, dipole,    &
                                 lowdin, mulliken, print_coeffs, number_restr, &
                                 Dbug, steep, Force_cut, Energy_cut, charge,   &
                                 minimzation_steep, n_min_steeps, n_points,    &
                                 lineal_search, timers, IGRID, IGRID2,         &
                                 use_libxc, ex_functional_id, ec_functional_id,&
                                 gpu_level, becke, PBE0
   use tbdft_data         , only: tbdft_calc, MTB, alfaTB, betaTB, gammaTB,    &
                                  start_tdtb, end_tdtb,n_biasTB,               &
                                  driving_rateTB, TB_q_tot, TB_charge_ref,     &
                                  TB_q_told
   use ECP_mod           , only: ecpmode, ecptypes, tipeECP, ZlistECP,         &
                                 verbose_ECP, cutECP, local_nonlocal,          &
                                 ecp_debug, FOCK_ECP_read, FOCK_ECP_write,     &
                                 ecp_full_range_int, Fulltimer_ECP, cut2_0,    &
                                 cut3_0
   use ehrendata         , only: ndyn_steps, edyn_steps, nullify_forces,       &
                                 wdip_nfreq, wdip_fname, rsti_loads,           &
                                 rsti_fname, rsto_saves, rsto_nfreq,           &
                                 rsto_fname, eefld_on, eefld_ampx, eefld_ampy, &
                                 eefld_ampz, eefld_timeamp, eefld_timegfh,     &
                                 eefld_timepos, eefld_timegih, eefld_wavelen
   use field_data        , only: field, a0, epsilon, Fx, Fy, Fz,               &
                                 field_iso_file, field_aniso_file,             &
                                 nfields_iso, nfields_aniso
   use fileio_data       , only: verbose, style, rst_dens, movie_nfreq,        &
                                 movie_name0
   use fockbias_data     , only: fockbias_is_active, fockbias_is_shaped,       &
                                 fockbias_timegrow , fockbias_timefall,        &
                                 fockbias_timeamp0 , fockbias_readfile
   use initial_guess_data, only: initial_guess
   use td_data           , only: tdrestart, writedens, td_rst_freq, tdstep,    &
                                 ntdstep, timedep, td_do_pop
   use trans_Data        , only: gaussian_convert
   use transport_data    , only: transport_calc, generate_rho0, nbias,         &
                                 save_charge_freq, driving_rate, Pop_Drive
   use ghost_atoms_data  , only: n_ghosts, ghost_atoms
   use basis_data        , only: norm, int_basis, rmax, rmaxs, basis_set,      &
                                 fitting_set
   use excited_data      , only: lresp, nstates, tolv, tole, fittExcited,      &
                                 libint_recalc, root, FCA, nfo, nfv, TSH,      &
                                 excited_forces
   use converger_data    , only: DIIS, ndiis, GOLD, told, Etold, good_cut,     &
                                 hybrid_converg, DIIS_bias, conver_method,     &
                                 level_shift, lvl_shift_cut, lvl_shift_en,     &
                                 Rho_LS, nMax, DIIS_start, BDIIS_start
   use dos_data          , only: dos_calc, pdos_calc, pdos_allb
   use dftd3_data        , only: dftd3
   use rhoint            , only: write_int_rho, w_rho_xmin, w_rho_ymin,        &
                                 w_rho_zmin, w_rho_xmax, w_rho_ymax,           &
                                 w_rho_zmax, w_rho_dx,  w_rho_dy, w_rho_dz,    &
                                 w_rho_rmin, w_rho_rmax, w_rho_dr,             &
                                 w_rho_dtheta, w_rho_dphi, write1Drho

   
   implicit none

!  Namelist definition
   namelist /lionml/ ndyn_steps, edyn_steps, nullify_forces, propagator,       &
                     wdip_nfreq, wdip_fname,                                   &
                     rsti_loads, rsti_fname, rsto_saves, rsto_nfreq,           &
                     rsto_fname,                                               &
                     eefld_on, eefld_ampx, eefld_ampy, eefld_ampz,             &
                     eefld_wavelen, eefld_timegih, eefld_timegfh,              &
                     eefld_timepos, eefld_timeamp,                             &
                     fockbias_is_active, fockbias_is_shaped, fockbias_readfile,&
                     fockbias_timegrow , fockbias_timefall , fockbias_timeamp0,&
                     use_libxc, ex_functional_id, ec_functional_id

   namelist /lio/ OPEN, NMAX, Nunp, VCINP, rmax, rmaxs, predcoef, writexyz,    &
                  Iexch, igrid, igrid2, initial_guess, natom, nsol, charge,    &
                  ! Convergence acceleration.
                  GOLD, told, Etold, good_cut, DIIS, ndiis, hybrid_converg,    &
                  diis_bias, conver_method, level_shift, lvl_shift_cut,        &
                  lvl_shift_en, DIIS_start, BDIIS_start,                       &
                  ! File Input/Output.
                  frestartin, style, frestart, fukui, dipole, lowdin, verbose, &
                  mulliken, writeforces, int_basis, fitting_set, basis_set,    &
                  restart_freq, print_coeffs, Dbug, timers, gaussian_convert,  &
                  rst_dens, becke,                                             &
                  ! DFT and TD-DFT Variables.
                  timedep, tdstep, ntdstep, propagator, NBCH, tdrestart,       &
                  writedens, td_rst_freq, td_do_pop,                           &
                  ! Field Variables
                  field, epsilon, a0, Fx, Fy, Fz, nfields_iso, nfields_aniso,  &
                  field_aniso_file, field_iso_file,                            &
                  ! Effective Core Potential Variables.
                  ecpmode, ecptypes, tipeECP, ZlistECP, cutECP,                &
                  local_nonlocal, ecp_debug, ecp_full_range_int, verbose_ECP,  &
                  FOCK_ECP_read, FOCK_ECP_write, Fulltimer_ECP, cut2_0, cut3_0,&
                  ! Geometry optimizations and restraints
                  steep, Force_cut, Energy_cut, minimzation_steep,             &
                  n_min_steeps,lineal_search,n_points, number_restr,           &
                  ! Variables for orbital printing.
                  cubegen_only, cube_res, cube_sel, cube_dens, cube_dens_file, &
                  cube_orb, cube_orb_file, cube_elec, cube_elec_file,          &
                  cube_sqrt_orb,                                               &
                  ! Variables for 1D density printing.
                  write_int_rho, w_rho_xmin, w_rho_ymin, w_rho_zmin,           &
                  w_rho_xmax, w_rho_ymax, w_rho_zmax, w_rho_dx,  w_rho_dy,     &
                  w_rho_dz, w_rho_rmin, w_rho_rmax, w_rho_dr, w_rho_dtheta,    &
                  w_rho_dphi, write1Drho,                                      &
                  ! Variables for GPU options.
                  little_cube_size, max_function_exponent, free_global_memory, &
                  min_points_per_cube, assign_all_functions, sphere_radius,    &
                  remove_zero_weights, energy_all_iterations, gpu_level,       &
                  ! Variables for Transport
                  transport_calc, generate_rho0, nbias,                        &
                  save_charge_freq, driving_rate, Pop_Drive,                   &
                  ! Variables for TBDFT
                  tbdft_calc, MTB, alfaTB, betaTB, gammaTB, start_tdtb,        &
                  end_tdtb,n_biasTB, driving_rateTB, TB_q_tot, TB_charge_ref,  &
                  TB_q_told,                                                   &
                  !Fockbias
                  fockbias_is_active, fockbias_is_shaped, fockbias_readfile,   &
                  fockbias_timegrow , fockbias_timefall , fockbias_timeamp0,   &
                  ! Libxc variables
                  use_libxc, ex_functional_id, ec_functional_id,               &
                  ! Variables for Ghost atoms:
                  n_ghosts, ghost_atoms,                                       &
                  ! Variables for Linear Response
                  lresp, nstates, tolv, tole, fittExcited, libint_recalc, root,&
                  FCA, nfo, nfv, TSH, excited_forces,                          &
                  ! linear search for rho
                  Rho_LS,                                                      &
                  !DOS-PDOS calc
                  dos_calc, pdos_calc, pdos_allb,                              &
                  ! Movie setups
                  movie_nfreq, movie_name0,                                    &
                  ! Dispersion corrections.
                  dftd3,                                                       &
                  ! PBE0 functional
                  PBE0

   type lio_input_data
      ! COMMON
      double precision :: etold, gold, good_cut, rmax, rmaxs, told, DIIS_bias, &
                          lvl_shift_cut, lvl_shift_en, DIIS_start, bDIIS_start
      integer          :: charge, iexch, igrid, igrid2, initial_guess, natom,  &
                          ndiis, nmax, nsol, nunp, conver_method
      logical          :: diis, hybrid_converg, open, predcoef, vcinp, &
                          writexyz, level_shift
      ! FILE IO
      character*20     :: frestartin, frestart
      character*40     :: basis_set, fitting_set
      integer          :: restart_freq, timers, verbose, rst_dens
      logical          :: dbug, dipole, fukui, gaussian_convert, int_basis,   &
                          lowdin, mulliken, print_coeffs, style, writeforces, &
                          becke
      ! TD-DFT and FIELD
      character*20     :: field_aniso_file, field_iso_file
      double precision :: a0, epsilon, Fx, Fy, Fz, tdstep
      integer          :: NBCH, nfields_aniso, nfields_iso, ntdstep,           &
                          propagator, td_rst_freq, timedep, td_do_pop
      logical          :: tdrestart, writedens, field
      ! ECP
      character*30     :: tipeECP
      double precision :: cut2_0, cut3_0
      integer          :: ecptypes, local_nonlocal, verbose_ECP, ZlistECP(128)
      logical          :: cutECP, ecp_debug, ecp_full_range_int, ecpmode,      &
                          FOCK_ECP_read, FOCK_ECP_write, Fulltimer_ECP
      ! Minimizations and restraints
      double precision :: Force_cut, Energy_cut, minimzation_steep
      integer          :: n_min_steeps, n_points, number_restr
      logical          :: lineal_search, steep
      ! CUBEGEN options
      character*20     :: cube_dens_file, cube_elec_file, cube_orb_file
      integer          :: cube_res,cube_sel
      logical          :: cube_dens, cube_elec, cube_orb, cubegen_only,        &
                          cube_sqrt_orb
      ! GPU Options
      double precision :: free_global_memory, little_cube_size, sphere_radius
      integer          :: min_points_per_cube, max_function_exponent, gpu_level
      logical          :: assign_all_functions, energy_all_iterations,         &
                          remove_zero_weights
      ! Transport and TBDFT
      double precision :: alfaTB, betaTB, driving_rate, gammaTB, Vbias_TB,     &
                          driving_rateTB, TB_charge_ref, TB_q_told
      logical          :: gate_field, generate_rho0, transport_calc
      integer          :: tbdft_calc, end_bTB, end_tdtb, MTB, pop_drive,       &
                          save_charge_freq, start_tdtb, nbias, n_biasTB,       &
                          TB_q_tot
      ! Ehrenfest
      character*80     :: rsti_fname, rsto_fname, wdip_fname
      double precision :: eefld_ampx, eefld_ampy, eefld_ampz, eefld_timeamp,   &
                          eefld_timepos, eefld_wavelen
      integer          :: edyn_steps, ndyn_steps, rsto_nfreq, wdip_nfreq
      logical          :: eefld_on, eefld_timegih, eefld_timegfh,              &
                          nullify_forces, rsti_loads, rsto_saves
      ! Fock Bias Potential
      character*80     :: fockbias_readfile
      double precision :: fockbias_timeamp0, fockbias_timefall,fockbias_timegrow
      logical          :: fockbias_is_active, fockbias_is_shaped

      ! Libxc configuration
      integer          :: ex_functional_id, ec_functional_id
      logical          :: use_libxc
      ! Ghost atoms
      integer          :: n_ghosts, ghost_atoms(300)
      !DOS-PDOS
      logical          :: dos_calc, pdos_calc, pdos_allb
      ! DFTD3
      logical          :: dftd3
   end type lio_input_data
contains

subroutine get_namelist(lio_in)
   implicit none
   type(lio_input_data), intent(out) :: lio_in

   ! General
   lio_in%etold          = etold         ; lio_in%gold       = gold
   lio_in%good_cut       = good_cut      ; lio_in%rmax       = rmax
   lio_in%rmaxs          = rmaxs         ; lio_in%told       = told
   lio_in%charge         = charge        ; lio_in%iexch      = iexch
   lio_in%igrid          = igrid         ; lio_in%igrid2     = igrid2
   lio_in%initial_guess  = initial_guess ; lio_in%natom      = natom
   lio_in%ndiis          = ndiis         ; lio_in%nmax       = nmax
   lio_in%nsol           = nsol          ; lio_in%nunp       = nunp
   lio_in%diis           = diis          ; lio_in%open       = open
   lio_in%hybrid_converg = hybrid_converg; lio_in%diis_bias  = diis_bias
   lio_in%conver_method  = conver_method ; lio_in%predcoef   = predcoef
   lio_in%level_shift    = level_shift   ; lio_in%diis_start = diis_start
   lio_in%lvl_shift_en   = lvl_shift_en  ; lio_in%bdiis_start= bdiis_start
   lio_in%lvl_shift_cut  = lvl_shift_cut ;

   ! Fileio
   lio_in%vcinp            = vcinp           ; lio_in%writexyz    = writexyz
   lio_in%frestartin       = frestartin      ; lio_in%frestart    = frestart
   lio_in%basis_set        = basis_set       ; lio_in%fitting_set = fitting_set
   lio_in%restart_freq     = restart_freq    ; lio_in%timers      = timers
   lio_in%dbug             = dbug            ; lio_in%dipole      = dipole
   lio_in%gaussian_convert = gaussian_convert; lio_in%fukui       = fukui
   lio_in%int_basis        = int_basis       ; lio_in%lowdin      = lowdin
   lio_in%mulliken         = mulliken        ; lio_in%style       = style
   lio_in%print_coeffs     = print_coeffs    ; lio_in%writeforces = writeforces
   lio_in%verbose          = verbose         ; lio_in%rst_dens    = rst_dens
   lio_in%becke            = becke           ;

   ! TDDFT - Fields
   lio_in%field_aniso_file = field_aniso_file; lio_in%a0         = a0
   lio_in%field_iso_file   = field_iso_file  ; lio_in%epsilon    = epsilon
   lio_in%nfields_aniso    = nfields_aniso   ; lio_in%Fx         = Fx
   lio_in%nfields_iso      = nfields_iso     ; lio_in%Fy         = Fy
   lio_in%td_rst_freq      = td_rst_freq     ; lio_in%Fz         = Fz
   lio_in%tdstep           = tdstep          ; lio_in%NBCH       = NBCH
   lio_in%ntdstep          = ntdstep         ; lio_in%propagator = propagator
   lio_in%timedep          = timedep         ; lio_in%tdrestart  = tdrestart
   lio_in%writedens        = writedens       ; lio_in%field      = field
   lio_in%td_do_pop        = td_do_pop       ;
   ! ECP
   lio_in%ecp_full_range_int = ecp_full_range_int; lio_in%cut2_0    = cut2_0
   lio_in%verbose_ECP        = verbose_ECP       ; lio_in%cut3_0    = cut3_0
   lio_in%local_nonlocal     = local_nonlocal    ; lio_in%ecptypes  = ecptypes
   lio_in%ZlistECP           = ZlistECP          ; lio_in%cutECP    = cutECP
   lio_in%Fulltimer_ECP      = Fulltimer_ECP     ; lio_in%ecp_debug = ecp_debug
   lio_in%FOCK_ECP_write     = FOCK_ECP_write    ; lio_in%ecpmode   = ecpmode
   lio_in%FOCK_ECP_read      = FOCK_ECP_read     ; lio_in%tipeECP   = tipeECP
   ! Minimization - Restraints
   lio_in%lineal_search     = lineal_search    ; lio_in%Energy_cut = Energy_cut
   lio_in%minimzation_steep = minimzation_steep; lio_in%steep      = steep
   lio_in%n_min_steeps      = n_min_steeps     ; lio_in%n_points   = n_points
   lio_in%number_restr      = number_restr     ; lio_in%Force_cut  = Force_cut
   ! CubeGen
   lio_in%cube_dens_file = cube_dens_file ; lio_in%cube_res  = cube_res
   lio_in%cube_elec_file = cube_elec_file ; lio_in%cube_sel  = cube_sel
   lio_in%cube_orb_file  = cube_orb_file  ; lio_in%cube_dens = cube_dens
   lio_in%cubegen_only   = cubegen_only   ; lio_in%cube_elec = cube_elec
   lio_in%cube_sqrt_orb  = cube_sqrt_orb  ; lio_in%cube_orb  = cube_orb
   ! GPU Options
   lio_in%free_global_memory    = free_global_memory
   lio_in%little_cube_size      = little_cube_size
   lio_in%sphere_radius         = sphere_radius
   lio_in%min_points_per_cube   = min_points_per_cube
   lio_in%max_function_exponent = max_function_exponent
   lio_in%assign_all_functions  = assign_all_functions
   lio_in%energy_all_iterations = energy_all_iterations
   lio_in%remove_zero_weights   = remove_zero_weights
   lio_in%gpu_level             = gpu_level
   ! Transport and TBDFT
   lio_in%driving_rate     = driving_rate    ; lio_in%alfaTB    = alfaTB
   lio_in%tbdft_calc        = tbdft_calc     ; lio_in%betaTB    = betaTB
   lio_in%nbias            = nbias           ; lio_in%gammaTB   = gammaTB
   lio_in%generate_rho0    = generate_rho0   ; lio_in%TB_q_tot  = TB_q_tot
   lio_in%transport_calc   = transport_calc  ; lio_in%n_biasTB  = n_biasTB
   lio_in%end_tdtb         = end_tdtb        ; lio_in%pop_drive = pop_drive
   lio_in%save_charge_freq = save_charge_freq; lio_in%MTB       = MTB
   lio_in%start_tdtb       = start_tdtb      ; lio_in%TB_q_told = TB_q_told
   lio_in%TB_charge_ref    = TB_charge_ref
   lio_in%driving_rateTB   = driving_rateTB
   ! Ghost atoms
   lio_in%n_ghosts = n_ghosts ; lio_in%ghost_atoms = ghost_atoms

   ! Ehrenfest
   lio_in%rsti_fname = rsti_fname; lio_in%eefld_timeamp  = eefld_timeamp
   lio_in%rsto_fname = rsto_fname; lio_in%eefld_timepos  = eefld_timepos
   lio_in%wdip_fname = wdip_fname; lio_in%eefld_wavelen  = eefld_wavelen
   lio_in%eefld_ampx = eefld_ampx; lio_in%eefld_timegih  = eefld_timegih
   lio_in%eefld_ampy = eefld_ampy; lio_in%eefld_timegfh  = eefld_timegfh
   lio_in%eefld_ampz = eefld_ampz; lio_in%nullify_forces = nullify_forces
   lio_in%edyn_steps = edyn_steps; lio_in%ndyn_steps     = ndyn_steps
   lio_in%rsto_nfreq = rsto_nfreq; lio_in%rsto_saves     = rsto_saves
   lio_in%wdip_nfreq = wdip_nfreq; lio_in%eefld_on       = eefld_on
   lio_in%rsti_loads = rsti_loads;
   ! Fock Bias Potential
   lio_in%fockbias_readfile  = fockbias_readfile
   lio_in%fockbias_timeamp0  = fockbias_timeamp0
   lio_in%fockbias_timefall  = fockbias_timefall
   lio_in%fockbias_timegrow  = fockbias_timegrow
   lio_in%fockbias_is_active = fockbias_is_active
   lio_in%fockbias_is_shaped = fockbias_is_shaped
   ! DOS-PDOS calc
   lio_in%dos_calc = dos_calc
   lio_in%pdos_calc= pdos_calc
   lio_in%pdos_allb= pdos_allb
   
   ! Dispersion corrections
   lio_in%dftd3 = dftd3
   ! Libxc configuration
   !lio_in%ex_functional_id = ex_functional_id
   !lio_in%ec_functional_id = ec_functional_id
   !lio_in%use_libxc = use_libxc

   return
end subroutine get_namelist

end module lionml_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
