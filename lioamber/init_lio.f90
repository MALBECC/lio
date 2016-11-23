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
                             energy_freq, style, allnml, writeforces,          &
                             cube_elec, cube_elec_file
      use ECP_mod   , only : ecpmode, ecptypes, tipeECP, ZlistECP, cutECP,     &
                             local_nonlocal, ecp_debug, ecp_full_range_int,    &
                             verbose_ECP, Cnorm, FOCK_ECP_read, FOCK_ECP_write,&
                             Fulltimer_ECP, cut2_0, cut3_0

      implicit none

!     Names of files used for input and output.
      basis          = 'basis'       ; output             = 'output'      ;
      fmulliken      = 'mulliken'    ; fcoord             = 'qm.xyz'      ;

!     Theory level options.
      OPEN           = .false.       ; told               = 1.0D-6        ;
      NMAX           = 100           ; Etold              = 1.0d0         ;
      basis_set      = "DZVP"        ; hybrid_converg     = .false.       ;
      int_basis      = .false.       ; good_cut           = 1D-5          ;
      DIIS           = .true.        ; rmax               = 16            ;
      ndiis          = 30            ; rmaxs              = 5             ;
      GOLD           = 10.           ; omit_bas           = .false.       ;
      fitting_set    = "DZVP Coulomb Fitting" ;

!     Effective Core Potential options.
      ecpmode        = .false.       ; cut2_0             = 15.d0         ;
      ecptypes       = 0             ; cut3_0             = 12.d0         ;
      tipeECP        = 'NOT-DEFINED' ; verbose_ECP        = 0             ;
      ZlistECP       = 0             ; ecp_debug          = .false.       ;
      FOCK_ECP_read  = .false.       ; Fulltimer_ECP      = .false.       ;
      FOCK_ECP_write = .false.       ; local_nonlocal     = 0             ;
      cutECP         = .true.        ; ecp_full_range_int = .false.       ;

!     TD-DFT options.
      timedep        = 0             ; Fx                 = 0.05          ;
      propagator     = 1             ; Fy                 = 0.05          ;
      tdstep         = 2.D-3         ; Fz                 = 0.05          ;
      ntdstep        = 1             ; tdrestart          = .false.       ;
      NBCH           = 10            ; exter              = .false.       ;
      field          = .false.       ;

!     Write options and Restart options.
      verbose        = .false.       ; writexyz           = .true.        ;
      writedens      = .false.       ; frestart           ='restart.out'  ;
      VCINP          = .false.       ; frestartin         = 'restart.in'  ;
      restart_freq   = 1             ; writeforces        = .false.       ;

!     Cube, grid and other options.
      predcoef       = .false.       ; cubegen_only       = .false.       ;
      idip           = 0             ; cube_res           = 40            ;
      intsoldouble   = .true.        ; cube_dens          = .false.       ;
      dgtrig         = 100.          ; cube_orb           = .false.       ;
      Iexch          = 9             ; cube_sel           = 0             ;
      integ          = .true.        ; cube_orb_file      = "orb.cube"    ;
      DENS           = .true.        ; cube_dens_file     = 'dens.cube'   ;
      IGRID          = 2             ; cube_elec          = .false.       ;
      IGRID2         = 2             ; cube_elec_file     = 'field.cube'  ;
      a0             = 1000.0        ; style              = .true.        ;
      epsilon        = 1.D0          ; allnml             = .true.        ;
      NUNP           = 0             ; energy_freq        = 1             ;

      return
      end subroutine lio_defaults
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% LIO_INIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Subroutine lio_defaults gives default values to LIO runtime options.         !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      subroutine init_lio_common(natomin, Izin, nclatom, charge, callfrom)

      use garcha_mod, only : idip, nunp, X, XX, RMM, d, c, a, Nuc, ncont, cx,  &
                             ax, Nucx, ncontx, cd, ad, Nucd, ncontd, indexii,  &
                             indexiid, r, v, rqm, Em, Rm, pc, nnat, af, B, Iz, &
                             natom, nco, ng0, ngd0, ngrid, nl, norbit, ntatom, &
                             allnml, style                             
      use ECP_mod,    only : Cnorm, ecpmode

      implicit none
      integer , intent(in) :: charge, nclatom, natomin, Izin(natomin), callfrom
      integer              :: i, ng2, ng3, ngdDyn, ngDyn, nqnuc, ierr, ios
!      integer              :: i, ng2, ng3, ngdnu, ngnu, ngdDyn, ngDyn, nqnuc,  &
!                              ierr, ios

!     call g2g_timer_start('lio_init')

! Some important values are: 
! Ngrid may be set to 0  in the case of Numerical Integration.                 !
! ngDyn  : n° of atoms times the n° of basis functions.                        !
! norbit : n° of molecular orbitals involved.                                  !
! ngdDyn : n° of atoms times the n° of auxiliary functions.                    !
! Ngrid  : n° of grid points (LS-SCF part).                                    !
! For large systems, ng2 may result in <0 due to overflow.                     !
      if (callfrom.eq.1) then
          natom  = natomin          ;  Iz = Izin  ;
          ntatom = natom + nclatom  ;
          allocate(r(ntatom,3), rqm(natom,3), pc(ntatom))
      endif
 


      call DIMdrive(ngDyn,ngdDyn) !found dimention for impoortant arrays

! falta cambiar ngnu y ngdnu
  
!      ngnu   = natom*ng0  ;  ngdnu  = natom*ngd0       ;
!      ngDyn  = ngnu       ;  ngdDyn = ngdnu            ;

      ng2    = 5*ngDyn*(ngDyn+1)/2 + 3*ngdDyn*(ngdDyn+1)/2 + ngDyn +           &
               ngDyn*norbit + Ngrid
      ng3    = 4*ngDyn

      allocate(X(ngDyn,ng3)  , XX(ngdDyn,ngdDyn) , RMM(ng2)    , d(natom, natom), &
               c(ngDyn,nl)   , a(ngDyn,nl)       , Nuc(ngDyn)  , ncont(ngDyn)   , &
               cx(ngdDyn,nl) , ax(ngdDyn,nl)     , Nucx(ngdDyn), ncontx(ngdDyn) , &
               cd(ngdDyn,nl) , ad(ngdDyn,nl)     , Nucd(ngdDyn), ncontd(ngdDyn) , &
               indexii(ngDyn), indexiid(ngdDyn)  , v(ntatom,3) , Em(ntatom)     , &
               Rm(ntatom)    , af(ngdDyn)         , nnat(ntatom), B(ngdDyn,3))

      ! Cnorm contains normalized coefficients of basis.
      ! Differentiate C for x^2,y^2,z^2 and  xy,xz,yx (3^0.5 factor)
      if (ecpmode) allocate (Cnorm(ngDyn,nl)) 

      call g2g_init()
      nqnuc = 0
      do i = 1, natom
          nqnuc = nqnuc + Iz(i)
      enddo
      nco = ((nqnuc - charge) - Nunp)/2

      ! Header for the file containing dipole moments.
      if (idip.eq.1) call write_dip_header(69)

!     Prints LIO logo to output and options chosen for the run. 
      if (style) call LIO_LOGO()
      if (style) call NEW_WRITE_NML(charge)
      call drive(ng2, ngDyn, ngdDyn)

!     call g2g_timer_stop('lio_init')
      return 
      end subroutine init_lio_common
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% INIT_LIO_GROMACS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Subroutine init_lio_gromacs performs Lio initialization when called from     !
! Gromacs software package, in order to conduct a hybrid QM/MM calculation.    !
! See also SCF_gro.f, which deals with the SCF cycle within the MD.            !
! In order to avoid compatibility problems due to a FORTRAN/C++ interface, LIO !
! options are read from a file named "lio.in" in the current workspace.        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      subroutine init_lio_gromacs(natomin, Izin, nclatom, chargein)

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
                             cube_elec, writeforces
      use ECP_mod   , only : ecpmode, ecptypes, tipeECP, ZlistECP, cutECP,     &
                             local_nonlocal, ecp_debug, ecp_full_range_int,    &
                             verbose_ECP, Cnorm, FOCK_ECP_read, FOCK_ECP_write,&
                             Fulltimer_ECP, cut2_0, cut3_0

      implicit none
      integer,  intent(in) :: chargein, nclatom, natomin, Izin(natomin)
      integer              :: ierr, ios 
      logical              :: file_exists
      character(len=20)    :: input_file

!                    Common LIO variables.
      namelist /lio/ OPEN, NMAX, Nunp, VCINP, frestartin, GOLD, told, Etold,   &
                     rmax, rmaxs, predcoef, idip, writexyz, intsoldouble, DIIS,&
                     ndiis, dgtrig, Iexch, integ, dens, igrid, igrid2,         &
                     hybrid_converg, good_cut, style, allnml, frestart,        &
!                    TD-DFT Variables.
                     timedep, tdstep, ntdstep, propagator, NBCH, field, a0,    &
                     epsilon, exter, Fx, Fy, Fz, tdrestart, writedens,         &
                     writeforces, basis_set, fitting_set, int_basis,           &
!                    Effective Core Potential Variables.
                     ecpmode, ecptypes, tipeECP,ZlistECP, cutECP,              &
                     local_nonlocal, ecp_debug,ecp_full_range_int,verbose_ECP, &
                     verbose, FOCK_ECP_read, FOCK_ECP_write, Fulltimer_ECP,    &
                     cut2_0, cut3_0,                                           &
!                    Cube variables.
                     cubegen_only, cube_res, cube_dens, cube_dens_file,        &
                     cube_orb, cube_sel, cube_orb_file, cube_elec,             &
                     cube_elec_file

!     Gives default values to variables.
      call lio_defaults()

!     Checks if input file exists and writes data to namelist variables.
      input_file = 'lio.in'
      inquire(file = input_file, exist = file_exists)
      if(file_exists) then
          open(unit = 100, file = input_file, iostat = ios)
          read(100, nml = lio, iostat = ierr)
          if(ierr.gt.0) stop 'Input error in LIO namelist.'
      else
          write(*,*) 'Inputfile lio.in not found. Using defaults.'
      endif

!     Initializes LIO. The last argument indicates LIO is not being used alone.
      call init_lio_common(natomin, Izin, nclatom, chargein, 1)

      return
      end subroutine init_lio_gromacs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% INIT_LIO_AMBER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Subroutine init_lio_amber performs Lio initialization when called from AMBER !
! software package, in order to conduct a hybrid QM/MM calculation. See also   !
! SCF_in.f, which deals with the SCF cycle within the MD.                      !
! AMBER directly passes options to LIO, but since the interface has not been   !
! officialy updated on the AMBER side, only some variables are received.       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
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

!     Gives default values to variables.       
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

!     Initializes LIO. The last argument indicates LIO is not being used alone.
      call init_lio_common(natomin, Izin, nclatom, charge, 1) 
      end subroutine init_lio_amber
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
