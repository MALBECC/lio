!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module lionml_data

   use garcha_mod        , only: natom, nsol, basis, output, fmulliken, fcoord,&
                                 OPEN, NMAX, basis_set, fitting_set, int_basis,&
                                 DIIS, ndiis, GOLD, told, Etold, good_cut,     &
                                 hybrid_converg, rmax, rmaxs, omit_bas, NBCH,  &
                                 propagator, verbose,VCINP, restart_freq,      &
                                 writexyz, dgtrig, Iexch, integ, frestartin,   &
                                 frestart, predcoef, idip, intsoldouble,       &
                                 cubegen_only, cube_res, cube_dens, cube_orb,  &
                                 DENS, cube_sel, cube_orb_file, cube_dens_file,&
                                 cube_elec, cube_elec_file, energy_freq, NUNP, &
                                 style, allnml, writeforces, cube_sqrt_orb,    &
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
   use field_data        , only: field, a0, epsilon, Fx, Fy, Fz,               &
                                 field_iso_file, field_aniso_file,             &
                                 nfields_iso, nfields_aniso
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

!  Run type information
!------------------------------------------------------------------------------!
   integer :: ndyn_steps = 0  ! Number of total nuclear dynamic steps
   integer :: edyn_steps = 0  ! Number of total electronic dynamic steps PER
                              ! nuclear dynamic step.
!
!  ndyn == 0 & edyn == 0   =>   Single point
!  ndyn /= 0 & edyn == 0   =>   BO atomistic dynamic (aux mm)
!  ndyn == 0 & edyn /= 0   =>   TD electron dynamic
!  ndyn /= 0 & edyn /= 0   =>   Ehrenfest dynamic (aux mm)
!
   logical :: nullify_forces = .false.
!  Output information
!------------------------------------------------------------------------------!
!  TODO: set option so that data is defaulted into one output and you
!        need to activate different files.
   integer           :: verbose_level = 0
   integer           :: wdip_nfreq = 0
   character(len=80) :: wdip_fname = "liorun.dip"
!
!
!  Restarting information
!------------------------------------------------------------------------------!
!
!  If (rsti_loads), the program will load the restart from (rsti_fname).
!
!  If (rsto_saves), the program will save the restart information overwriting
!  the file (rsto_fname) every (rsto_nfreq) steps and in the last step.
!
   logical           :: rsti_loads = .false.
   character(len=80) :: rsti_fname = "liorun.rsti"

   logical           :: rsto_saves = .false.
   integer           :: rsto_nfreq = 0
   character(len=80) :: rsto_fname = "liorun.rsto"
!
!  External Electrical Field
!------------------------------------------------------------------------------!
!     If (eefld_on), an external field will be applied to the system. The
!  amplitude in each direction is given by the (eefld_amp) variables. It
!  can have an oscilating time modulation of a specific (eefld_wavelen) and
!  also a gaussian envelop centered in (eefld_timepos), with width given by
!  (eefld_timeamp). Both (eefld_timegih) and (eefld_timegfh) must be true for
!  a full gaussian, activating the modulation before and after the center
!  respectively.
!
   logical :: eefld_on   = .false.
   real*8  :: eefld_ampx = 0.0d0 ! in au
   real*8  :: eefld_ampy = 0.0d0 ! in au
   real*8  :: eefld_ampz = 0.0d0 ! in au

   logical :: eefld_timegih = .false. ! time gaussian initial half
   logical :: eefld_timegfh = .false. ! time gaussian final half
   real*8  :: eefld_timepos =  1.0d0  ! in ps (currently fs!)
   real*8  :: eefld_timeamp =  0.2d0  ! in ps (currently fs!)
   real*8  :: eefld_wavelen =  0.0d0  ! in nm
!
!  Namelist definition
!------------------------------------------------------------------------------!
   namelist /lionml/ ndyn_steps, edyn_steps, nullify_forces, propagator,       &
                     verbose_level, wdip_nfreq, wdip_fname,                    &
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
                  initial_guess,                                               &
                  ! File Input/Output.
                  frestartin, style, allnml, frestart, fukui, dipole, lowdin,  &
                  mulliken, writeforces, int_basis, fitting_set, basis_set,    &
                  restart_freq, print_coeffs,                                  &
                  ! DFT and TD-DFT Variables.
                  timedep, tdstep, ntdstep, propagator, NBCH, tdrestart,       &
                  writedens, td_rst_freq,                                      &
                  ! Field Variables
                  field, epsilon, a0, Fx, Fy, Fz, nfields_iso, nfields_aniso,  &
                  field_aniso_file, field_iso_file,                            &
                  ! Effective Core Potential Variables.
                  ecpmode, ecptypes, tipeECP, ZlistECP, cutECP, ecp_debug,     &
                  local_nonlocal, ecp_debug, ecp_full_range_int, verbose_ECP,  &
                  verbose, FOCK_ECP_read, FOCK_ECP_write, Fulltimer_ECP,       &
                  cut2_0, cut3_0,                                              &
                  ! Distance Restrain
                  number_restr,                                                &
                  ! Debug variables
                  Dbug, timers,                                                &
                  ! Geometry optimizations
                  steep, Force_cut, Energy_cut, minimzation_steep,             &
                  n_min_steeps,lineal_search,n_points,                         &
                  ! Variables for orbital printing.
                  cubegen_only, cube_res, cube_sel, cube_dens, cube_dens_file, &
                  cube_orb, cube_orb_file, cube_elec, cube_elec_file,          &
                  cube_sqrt_orb,                                               &
                  ! Variables for GPU options.
                  little_cube_size, max_function_exponent, free_global_memory, &
                  min_points_per_cube, assign_all_functions, sphere_radius,    &
                  remove_zero_weights, energy_all_iterations,                  &
                  ! Variables when LIO is used alone.
                  natom, nsol, charge,                                         &
                  ! Variables for Transport
                  transport_calc, generate_rho0, gate_field,                   &
                  save_charge_freq, driving_rate, Pop_Drive,                   &
                  ! Variables for DFTB
                  dftb_calc, MTB, alfaTB, betaTB, gammaTB, Vbias_TB, end_bTB,  &
                  start_tdtb, end_tdtb, TBsave, TBload,                        &
                  ! Variables for translation
                  gaussian_convert

end module lionml_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
