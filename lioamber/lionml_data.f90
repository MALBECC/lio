!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module lionml_data

   use garcha_mod        , only: natom, nsol, basis, output, fmulliken, fcoord,&
                                 OPEN, NMAX, basis_set, fitting_set, int_basis,&
                                 DIIS, ndiis, GOLD, told, Etold, good_cut,     &
                                 hybrid_converg, rmax, rmaxs, omit_bas, NBCH,  &
                                 propagator, VCINP, restart_freq,              &
                                 writexyz, dgtrig, Iexch, integ, frestartin,   &
                                 frestart, predcoef, idip, intsoldouble,       &
                                 cubegen_only, cube_res, cube_dens, cube_orb,  &
                                 DENS, cube_sel, cube_orb_file, cube_dens_file,&
                                 cube_elec, cube_elec_file, energy_freq, NUNP, &
                                 writeforces, cube_sqrt_orb,    &
                                 fukui, little_cube_size, min_points_per_cube, &
                                 max_function_exponent, assign_all_functions,  &
                                 remove_zero_weights, energy_all_iterations,   &
                                 free_global_memory, sphere_radius, dipole,    &
                                 lowdin, mulliken, print_coeffs, number_restr, &
                                 Dbug, steep, Force_cut, Energy_cut, charge,   &
                                 minimzation_steep, n_min_steeps, n_points,    &
                                 lineal_search, timers, IGRID, IGRID2
   use dftb_data         , only: dftb_calc, MTB, alfaTB, betaTB, gammaTB,      &
                                 Vbias_TB, end_bTB, start_tdtb, end_tdtb,      &
                                 TBsave, TBload
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
   use fileio_data       , only: verbose, style
   use fockbias_data     , only: fockbias_is_active, fockbias_is_shaped,       &
                                 fockbias_timegrow , fockbias_timefall,        &
                                 fockbias_timeamp0 , fockbias_readfile
   use initial_guess_data, only: initial_guess
   use td_data           , only: tdrestart, writedens, td_rst_freq, tdstep,    &
                                 ntdstep, timedep
   use trans_Data        , only: gaussian_convert
   use transport_data    , only: transport_calc, generate_rho0, gate_field,    &
                                 save_charge_freq, driving_rate, Pop_Drive
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
                     fockbias_timegrow , fockbias_timefall , fockbias_timeamp0

   namelist /lio/ OPEN, NMAX, Nunp, VCINP, GOLD, told, Etold, rmax, rmaxs,     &
                  predcoef, idip, writexyz, intsoldouble, DIIS, ndiis, dgtrig, &
                  Iexch, integ, dens, igrid, igrid2, good_cut, hybrid_converg, &
                  initial_guess, natom, nsol, charge,                          &
                  ! File Input/Output.
                  frestartin, style, frestart, fukui, dipole, lowdin,          &
                  mulliken, writeforces, int_basis, fitting_set, basis_set,    &
                  restart_freq, print_coeffs, Dbug, timers, gaussian_convert,  &
                  ! DFT and TD-DFT Variables.
                  timedep, tdstep, ntdstep, propagator, NBCH, tdrestart,       &
                  writedens, td_rst_freq,                                      &
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
                  ! Variables for GPU options.
                  little_cube_size, max_function_exponent, free_global_memory, &
                  min_points_per_cube, assign_all_functions, sphere_radius,    &
                  remove_zero_weights, energy_all_iterations,                  &
                  ! Variables for Transport
                  transport_calc, generate_rho0, gate_field,                   &
                  save_charge_freq, driving_rate, Pop_Drive,                   &
                  ! Variables for DFTB
                  dftb_calc, MTB, alfaTB, betaTB, gammaTB, Vbias_TB, end_bTB,  &
                  start_tdtb, end_tdtb, TBsave, TBload

   type lio_input_data
      ! COMMON
      double precision :: etold, gold, good_cut, rmax, rmaxs, told
      integer          :: charge, iexch, igrid, igrid2, initial_guess, natom,  &
                          ndiis, nmax, nsol, nunp
      logical          :: dens, diis, hybrid_converg, integ, intsoldouble,     &
                          open, predcoef, vcinp, writexyz
      ! FILE IO
      character*20     :: frestartin, frestart
      character*40     :: basis_set, fitting_set
      integer          :: restart_freq, timers
      logical          :: dbug, dipole, fukui, gaussian_convert, int_basis,    &
                          lowdin, mulliken, print_coeffs, style, writeforces
      ! TD-DFT and FIELD
      character*20     :: field_aniso_file, field_iso_file
      double precision :: a0, epsilon, Fx, Fy, Fz, tdstep
      integer          :: NBCH, nfields_aniso, nfields_iso, ntdstep,           &
                          propagator, td_rst_freq, timedep
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
      integer          :: min_points_per_cube, max_function_exponent
      logical          :: assign_all_functions, energy_all_iterations,         &
                          remove_zero_weights
      ! Transport and DFTB
      double precision :: alfaTB, betaTB, driving_rate, gammaTB, Vbias_TB
      logical          :: dftb_calc, gate_field, generate_rho0, transport_calc
      integer          :: end_bTB, end_tdtb, MTB, pop_drive, save_charge_freq, &
                          start_tdtb, TBload, TBsave
   end type lio_input_data

   type lionml_input_data

   end type lionml_input_data
contains
end module lionml_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
