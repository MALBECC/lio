      subroutine init_lio_amber(natomin, Izin, nclatom, charge, basis_i        &
                 , output_i, fcoord_i, fmulliken_i, frestart_i, frestartin_i   &
                 , verbose_i, OPEN_i, NMAX_i, NUNP_i, VCINP_i, GOLD_i, told_i  &
                 , rmax_i, rmaxs_i, predcoef_i, idip_i, writexyz_i             & 
                 , intsoldouble_i, DIIS_i, ndiis_i, dgtrig_i, Iexch_i, integ_i &
                 , DENS_i , IGRID_i, IGRID2_i , timedep_i , tdstep_i           &
                 , ntdstep_i, field_i, exter_i, a0_i, epsilon_i, Fx_i, Fy_i    &
                 , Fz_i, NBCH_i, propagator_i, writedens_i, tdrestart_i        &
#ifdef MOD_AMBER
                 , basis_set_i, fitting_set_i, int_basis_i, cubegen_only_i     &
                 , cuberes_i, cubedens_i, cubedensfile_i, cubeorb_i, cubesel_i &
                 , cubeorbfile_i, restart_freq_i, energy_freq_i)
#else
                 )
#endif

      use garcha_mod, only : basis, output, fmulliken, fcoord, OPEN, NMAX,     &
                             basis_set, fitting_set, int_basis, DIIS, ndiis,   &
                             GOLD, told, Etold, hybrid_converg, good_cut,      &
                             rmax, rmaxs, omit_bas, timedep, propagator,       &
                             tdstep, ntdstep, NBCH, Fx, Fy, Fz, field,         &
                             tdrestart, exter, verbose, writedens, VCINP,      &
                             restart_freq, writexyz, frestartin,               &
                             frestart, predcoef, idip, intsoldouble, dgtrig,   &
                             Iexch, integ, DENS, IGRID, IGRID2, a0, epsilon,   &
                             cubegen_only, cube_res, cube_dens, cube_orb,      &
                             cube_sel, cube_orb_file, cube_dens_file, NUNP,    &
                             energy_freq, style, allnml, cube_elec_file,       &
                             cube_elec
      use ECP_mod   , only : ecpmode, ecptypes, tipeECP, ZlistECP, cutECP,     &
                             local_nonlocal, ecp_debug, ecp_full_range_int,    &
                             verbose_ECP, Cnorm, FOCK_ECP_read, FOCK_ECP_write,&
                             Fulltimer_ECP, cut2_0, cut3_0

      implicit none
      integer , intent(in) :: charge, nclatom, natomin, Izin(natomin)
#ifdef MOD_AMBER
      character(len=40) :: basis_set_i, fitting_set_i
      character(len=20) :: cubedensfile_i,cubeorbfile_i
      logical           :: int_basis_i, cubegen_only_i, cubedens_i, cubeorb_i
      integer           :: cuberes_i, cubesel_i, restart_freq_i, energy_freq_i
#endif
      character(len=20) :: basis_i, output_i, fcoord_i, fmulliken_i,           &
                           frestart_i, frestartin_i
      logical           :: verbose_i, OPEN_i, VCINP_i, predcoef_i, writexyz_i, &
                           intsoldouble_i, DIIS_i, integ_i, DENS_i, field_i,   &
                           exter_i, writedens_i, tdrestart_i
      integer           :: NMAX_i, NUNP_i, idip_i, ndiis_i, Iexch_i, IGRID_i,  &
                           IGRID2_i, timedep_i, ntdstep_i, NBCH_i, propagator_i
      real*8            :: GOLD_i, told_i, rmax_i, rmaxs_i, dgtrig_i, tdstep_i,&
                           a0_i, epsilon_i, Fx_i, Fy_i, Fz_i
       
      call lio_defaults()  

#ifdef MOD_AMBER
      basis_set      = basis_set_i    ; fitting_set   = fitting_set_i  ;
      int_basis      = int_basis_i    ; cubegen_only  = cubegen_only_i ;
      cube_res       = cuberes_i      ; cube_dens     = cubedens_i     ;
      cube_dens_file = cubedensfile_i ; cube_orb      = cubeorb_i      ;
      cube_sel       = cubesel_i      ; cube_orb_file = cubeorbfile_i  ;
      restart_freq   = restart_freq_i ; energy_freq   = energy_freq_i  ;
#endif
      basis          = basis_i        ; output        = output_i       ;
      fcoord         = fcoord_i       ; fmulliken     = fmulliken_i    ;
      frestart       = frestart_i     ; frestartin    = frestartin_i   ;
      verbose        = verbose_i      ; OPEN          = OPEN_i         ;
      NMAX           = NMAX_i         ; NUNP          = NUNP_i         ;
      VCINP          = VCINP_i        ; GOLD          = GOLD_i         ;
      told           = told_i         ; rmax          = rmax_i         ;
      rmaxs          = rmaxs_i        ; predcoef      = predcoef_i     ;
      idip           = idip_i         ; writexyz      = writexyz_i     ;
      intsoldouble   = intsoldouble_i ; DIIS          = DIIS_i         ;
      ndiis          = ndiis_i        ; dgtrig        = dgtrig_i       ;
      Iexch          = Iexch_i        ; integ         = integ_i        ;
      DENS           = DENS_i         ; IGRID         = IGRID_i        ;
      IGRID2         = IGRID2_i       ; timedep       = timedep_i      ;
      field          = field_i        ; exter         = exter_i        ;
      tdstep         = tdstep_i       ; ntdstep       = ntdstep_i      ;
      a0             = a0_i           ; epsilon       = epsilon_i      ;
      Fx             = Fx_i           ; Fy            = Fy_i           ;
      Fz             = Fz_i           ; NBCH          = NBCH_i         ;
      propagator     = propagator_i   ; writedens     = writedens_i    ;
      tdrestart      = tdrestart_i

      call lio_init() 
      end
