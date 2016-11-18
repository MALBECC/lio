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
      
      call lio_init(natomin, Izin, nclatom, chargein)

      end
      
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

