!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% INIT_LIO_GROMACS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! File init_gromacs.f contains init_lio_gromacs, which performs Lio            !
! initialization when called from Gromacs software package in order to perform !
! a hybrid QM/MM calculation. See also SCF_gro.f, which deals with the SCF     !
! cycle within the molecular dynamics.                                         !
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
                             energy_freq, style, allnml, X, XX, RMM, d, c, a,  &
                             Nuc, ncont, cx, ax, Nucx, ncontx, cd, ad, Nucd,   &
                             ncontd, indexii, indexiid, r, v, rqm, Em, Rm, pc, &
                             nnat, af, B, Iz, natom, nco, ng0, ngd0, ngrid, nl,&
                             norbit, ntatom, cube_elec, cube_elec_file
      use ECP_mod   , only : ecpmode, ecptypes, tipeECP, ZlistECP, cutECP,     &
                             local_nonlocal, ecp_debug, ecp_full_range_int,    &
                             verbose_ECP, Cnorm, FOCK_ECP_read, FOCK_ECP_write,&
                             Fulltimer_ECP, cut2_0, cut3_0

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% Parameter Definition and Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      implicit none
      integer , intent(in) :: chargein, nclatom, natomin, Izin(natomin)
      integer              :: i, ng2, ng3, ngdnu, ngnu, ngdDyn, ngDyn, nqnuc,  &
                              ierr, ios, charge 
      logical              :: writeforces, file_exists
      character(len=20)    :: input_file

!                    Common LIO variables.
      namelist /lio/ OPEN, NMAX, Nunp, VCINP, frestartin, GOLD, told, Etold,   &
                     rmax, rmaxs, predcoef, idip, writexyz, intsoldouble, DIIS,&
                     ndiis, dgtrig, Iexch, integ, dens, igrid, igrid2,         &
                     hybrid_converg, good_cut, style, allnml,                  &
!                    TD-DFT Variables.
                     timedep, tdstep, ntdstep, propagator, NBCH, field, a0,    &
                     epsilon, exter, Fx, Fy, Fz, tdrestart, writedens,         &
                     writeforces, basis_set, fitting_set, int_basis,           &
!                    Effective Core Potential Variables.
                     ecpmode, ecptypes, tipeECP,ZlistECP, cutECP,              &
                     local_nonlocal, ecp_debug,ecp_full_range_int,verbose_ECP, &
                     verbose, FOCK_ECP_read, FOCK_ECP_write, Fulltimer_ECP,    &
                     cut2_0, cut3_0,                                           &
!                    Cube variables
                     cubegen_only, cube_res, cube_dens, cube_dens_file,        &
                     cube_orb, cube_sel, cube_orb_file, cube_elec,             &
                     cube_elec_file

!     Names of files used for input and output.
      basis     = 'basis'    ; output = 'output' ;
      fmulliken = 'mulliken' ; fcoord = 'qm.xyz' ;

!     Theory level options.
      OPEN      = .false. ; ndiis = 30     ; good_cut = 1D-5    ;
      NMAX      = 100     ; GOLD  = 10.    ; rmax     = 16      ;
      basis_set = "DZVP"  ; told  = 1.0D-6 ; rmaxs    = 5       ;
      fitting_set = "DZVP Coulomb Fitting" ; omit_bas = .false. ;
      int_basis = .false. ; Etold = 1.0d0  ;
      DIIS      = .true.  ; hybrid_converg = .false.            ; 

!     Effective Core Potential options.
      ecpmode        = .false.       ; cut2_0             = 15.d0   ;
      ecptypes       = 0             ; cut3_0             = 12.d0   ;
      tipeECP        = 'NOT-DEFINED' ; verbose_ECP        = 0       ;
      ZlistECP       = 0             ; ecp_debug          = .false. ;
      FOCK_ECP_read  = .false.       ; Fulltimer_ECP      = .false. ;
      FOCK_ECP_write = .false.       ; local_nonlocal     = 0       ;
      cutECP         = .true.        ; ecp_full_range_int = .false. ;

!     TD-DFT options.
      timedep    = 0     ; NBCH = 10   ; field     = .false. ;
      propagator = 1     ; Fx   = 0.05 ; tdrestart = .false. ;
      tdstep     = 2.D-3 ; Fy   = 0.05 ; exter     = .false. ;
      ntdstep    = 1     ; Fz   = 0.05 ;

!     Write options and Restart options.
      verbose   = .true.  ; writexyz    = .true.        ;
      writedens = .false. ; frestart    = 'restart.out' ;
      VCINP     = .false. ; frestartin  = 'restart.in'  ;
      restart_freq = 1    ; writeforces = .false.       ;

!     Cube, grid and other options.
      predcoef     = .false. ; cubegen_only   = .false.      ;
      idip         = 1       ; cube_res       = 40           ;
      intsoldouble = .true.  ; cube_dens      = .false.      ;
      dgtrig       = 100.    ; cube_orb       = .false.      ;
      Iexch        = 9       ; cube_sel       = 0            ;
      integ        = .true.  ; cube_orb_file  = "orb.cube"   ;
      DENS         = .true.  ; cube_dens_file = 'dens.cube'  ;
      IGRID        = 2       ; cube_elec      = .false.      ;
      IGRID2       = 2       ; cube_elec_file = 'field.cube' ;
      a0           = 1000.0  ; style          = .true.       ; 
      epsilon      = 1.D0    ; allnml         = .true.       ;
      NUNP         = 1       ; energy_freq    = 1            ;

!     Checks if input file exists and writes data to namelist variables.
      inquire(file = input_file, exist = file_exists)

      if(file_exists) then
          open(unit = 100, file = input_file, iostat = ios)
          read(100, nml = lio, iostat = ierr)
          if(ierr.gt.0) stop 'Input error in LIO namelist.'
      else
          write(*,*) 'Inputfile lio.in not found. Using defaults.'
      endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% Additional parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Natom is updated since it may change in every simulation step. Ngrid may be  !
! set to 0  in the case of Numerical Integration.                              !
! ngDyn  : n° of atoms times the n° of basis functions.                        !
! norbit : n° of molecular orbitals involved.                                  !
! ngdDyn : n° of atoms times the n° of auxiliary functions.                    !
! Ngrid  : n° of grid points (LS-SCF part).                                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      charge = chargein
      natom  = natomin         
      ngnu   = natom*ng0 
      ngdnu  = natom*ngd0  
      Iz     = Izin            
      ngDyn  = ngnu      
      ngdDyn = ngdnu      
      ntatom = natom + nclatom 
      ng2    = 5*ngDyn*(ngDyn+1)/2 + 3*ngdDyn*(ngdDyn+1)/2 + ngDyn +     &
               ngDyn*norbit + Ngrid 
      ng3    = 4*ngDyn         

      allocate(X(ngDyn,ng3) , XX(ngdDyn,ngdDyn), RMM(ng2)   , d(natom, natom), &
               c(ngnu,nl)   , a(ngnu,nl)       , Nuc(ngnu)  , ncont(ngnu)    , &
               cx(ngdnu,nl) , ax(ngdnu,nl)     , Nucx(ngdnu), ncontx(ngdnu)  , &
               cd(ngdnu,nl) , ad(ngdnu,nl)     , Nucd(ngdnu), ncontd(ngdnu)  , &
               indexii(ngnu), indexiid(ngdnu)  , r(ntatom,3), v(ntatom,3)    , &
               rqm(natom,3) , Em(ntatom)       , Rm(ntatom) , pc(ntatom)     , &
               nnat(ntatom) , af(natom*ngd0)   , B(natom*ngd0,3))

      call g2g_init()
      nqnuc = 0
      do i = 1, natom
          nqnuc = nqnuc + Iz(i)
      enddo
      nco = ((nqnuc - charge) - Nunp)/2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% Final Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      call LIO_LOGO()
      call NEW_WRITE_NML(charge)
      call drive(ng2, ngDyn, ngdDyn)
      write(*,*)     '--- LIO Initialization OK ---'
      end
      
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

