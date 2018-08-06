!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%  LIONML_SUBS.F90  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains procedures to handle both lio and lionml namelists. It    !
! includes the following subroutines:                                          !
! * lionml_read: Reads both lio and lionml namelists from input file.          !
! * lionml_write: Prints both lio and lionml namelists to standard output.     !
! * lionml_check: Performs consistency checks on the namelist keywords. (TO-DO)!
!                                                                              !
! In addition, the following subroutines are meant only accessed internally:   !
! * lionml_write_dull: Prints namelists in a straightforward manner.           !
! * lionml_write_style: Prints namelists in a fancy manner.                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module lionml_subs
   implicit none
contains

!  The parameters are self explanatory: file_unit is the unit of an already
!  opened file, and return stat is a status that lets the caller know wether
!  the namelist was correctly read/written/checked or if there was a problem
!  during these processes. For sake of simplification, the opening and closing
!  of input/output files must be handled externally.
!  TODO: Implement a (simple) logger and use it for the error messages.
subroutine lionml_check(extern_stat)
   implicit none
   integer, intent(out), optional :: extern_stat
   integer                        :: intern_stat

   intern_stat = 0
   if ( present(extern_stat) ) extern_stat = 0
   return
end subroutine lionml_check

subroutine lionml_read(file_unit, extern_stat )

   use lionml_data, only: lionml, lio
   use fileio_data, only: verbose
   implicit none
   integer, intent(in)            :: file_unit
   integer, intent(out), optional :: extern_stat
   integer                        :: intern_stat

   ! Old lio namelist.
   intern_stat = 0
   rewind( unit = file_unit, iostat = intern_stat )
   if ( intern_stat /= 0 ) then
      write(*,'(A)') &
         "Cannot rewind LIO input file. Using defaults for namelist lio."
      if (verbose .gt. 3) write(*,'(A,I4)') "iostat = ", intern_stat
      if ( present(extern_stat) ) extern_stat = 1
      return
   end if

   intern_stat = 0
   read( unit = file_unit, nml = lio, iostat = intern_stat )
   if ( intern_stat /= 0 ) then
      write(*,'(A)') &
         "Error found in lio namelist. Using defaults for namelist lio."
      if (verbose .gt. 3) write(*,'(A,I4)') "iostat = ", intern_stat
      if ( present(extern_stat) ) extern_stat = 2
   end if

   ! New lionml namelist.
   intern_stat = 0
   rewind( unit = file_unit, iostat = intern_stat )
   if ( intern_stat /= 0 ) then
      write(*,'(A)') &
         "Cannot rewind LIO input file. Using defaults for namelist lionml."
      if (verbose .gt. 3) write(*,'(A,I4)') "iostat = ", intern_stat
      if ( present(extern_stat) ) extern_stat = 1
      return
   end if

   intern_stat = 0
   read( unit = file_unit, nml = lionml, iostat = intern_stat )
   if ( intern_stat /= 0 ) then
      write(*,'(A)') &
         "Error found in lionml namelist. Using defaults for namelist lionml."
      if (verbose .gt. 3) write(*,'(A,I4)') "iostat = ", intern_stat
      if ( present(extern_stat) ) extern_stat = 2
      return
   end if

   if ( present(extern_stat) ) extern_stat = 0
   return
end subroutine lionml_read


subroutine lionml_write(extern_stat)
   use fileio_data, only: get_style
   implicit none
   integer, intent(out), optional :: extern_stat
   logical                        :: my_style

   call get_style(my_style)
   if (my_style) then
      call lionml_write_style()
   else
      call lionml_write_dull()
   endif

   return
end subroutine lionml_write

subroutine lionml_write_dull()
   use lionml_data, only: lio_input_data, get_namelist
   type(lio_input_data) :: inputs

   call get_namelist(inputs)
   if (inputs%verbose .lt. 1) return

   write(*,*)
   write(*,9000) "LIO input variables"
   write(*,9000) " ! -- General and theory level: -- !"
   write(*,8000) inputs%natom, inputs%nsol, inputs%charge, inputs%Nunp, &
                 inputs%open
   write(*,8001) inputs%nmax, inputs%int_basis, inputs%basis_set
   write(*,8002) inputs%fitting_set, inputs%diis, inputs%ndiis
   write(*,8003) inputs%hybrid_converg, inputs%gold, inputs%told
   write(*,8004) inputs%Etold, inputs%good_cut, inputs%Rmax
   write(*,8005) inputs%RmaxS, inputs%Iexch, inputs%Igrid, inputs%Igrid2
   write(*,8006) inputs%PredCoef, inputs%initial_guess, inputs%dbug
   write(*,8007) inputs%n_ghosts
   write(*,9000) " ! -- Input and output options: -- !"
   write(*,8020) inputs%verbose, inputs%style, inputs%timers, inputs%writexyz, &
                 inputs%WriteForces
   write(*,8021) inputs%dipole, inputs%mulliken, inputs%lowdin, inputs%fukui,  &
                 inputs%print_coeffs
   write(*,8022) inputs%vcinp, inputs%Frestartin, inputs%restart_freq
   write(*,8023) inputs%frestart, inputs%Tdrestart, inputs%writedens
   write(*,8024) inputs%td_rst_freq, inputs%gaussian_convert
   write(*,8025) inputs%rst_dens
   write(*,9000) " ! -- TD-DFT and external field: -- !"
   write(*,8040) inputs%timedep, inputs%ntdstep, inputs%tdstep,inputs%propagator
   write(*,8041) inputs%NBCH, inputs%field, inputs%a0, inputs%epsilon
   write(*,8042) inputs%Fx, inputs%Fy, inputs%Fz, inputs%nfields_iso
   write(*,8043) inputs%nfields_aniso, inputs%field_iso_file
   write(*,8044) inputs%field_aniso_file, inputs%td_do_pop
   write(*,9000) " ! -- Effective Core Potentials: -- !"
   write(*,8060) inputs%Ecpmode, inputs%Ecptypes, inputs%TipeECP
   write(*,8061) inputs%Fock_ECP_read, inputs%Fock_ECP_write, inputs%cutECP, &
                 inputs%cut2_0
   write(*,8062) inputs%cut3_0, inputs%Verbose_ECP, inputs%ECP_debug, &
                 inputs%fulltimer_ECP
   write(*,8063) inputs%local_nonlocal, inputs%ECP_full_range_int
   call write_Zlist_ECP_dull(inputs%ZlistECP, inputs%Ecptypes)
   write(*,9000) " ! -- Minimizations and restraints: -- !"
   write(*,8080) inputs%steep, inputs%minimzation_steep, inputs%Energy_cut
   write(*,8081) inputs%Force_cut, inputs%n_min_steeps, inputs%lineal_search
   write(*,8082) inputs%n_points, inputs%number_restr
   write(*,9000) " ! -- CubeGen: -- !"
   write(*,8100) inputs%Cubegen_only, inputs%Cube_Res, inputs%Cube_Sel, &
                 inputs%Cube_Dens
   write(*,8101) inputs%Cube_Orb, inputs%Cube_Elec, inputs%Cube_Dens_file
   write(*,8102) inputs%cube_orb_file, inputs%cube_dens_file
   write(*,9000) " ! -- GPU Options: -- !"
   write(*,8120) inputs%energy_all_iterations, inputs%assign_all_functions, &
                 inputs%remove_zero_weights
   write(*,8121) inputs%max_function_exponent, inputs%free_global_memory
   write(*,8122) inputs%sphere_radius, inputs%little_cube_size
   write(*,8123) inputs%min_points_per_cube
   write(*,9000) " ! -- Transport and DFTB: -- !"
   write(*,8140) inputs%transport_calc, inputs%generate_rho0, &
                 inputs%driving_rate
   write(*,8141) inputs%gate_field, inputs%pop_drive, inputs%save_charge_freq, &
                 inputs%dftb_calc
   write(*,8142) inputs%MTB, inputs%alfaTB, inputs%betaTB
   write(*,8143) inputs%gammaTB, inputs%Vbias_TB, inputs%start_tdtb
   write(*,8144) inputs%end_tdtb, inputs%end_bTB, inputs%TBload, inputs%TBsave
   write(*,8145) inputs%nbias
   write(*,9000) " ! -- Ehrenfest: -- !"
   write(*,8160) inputs%ndyn_steps, inputs%edyn_steps, inputs%nullify_forces
   write(*,8161) inputs%wdip_nfreq, inputs%wdip_fname, inputs%rsti_loads
   write(*,8162) inputs%rsto_saves, inputs%rsto_nfreq, inputs%rsti_fname
   write(*,8163) inputs%rsto_fname, inputs%eefld_on, inputs%eefld_ampx
   write(*,8164) inputs%eefld_ampy, inputs%eefld_ampz
   write(*,8165) inputs%eefld_timeamp, inputs%eefld_timepos
   write(*,8166) inputs%eefld_timegfh, inputs%eefld_timegih, &
                 inputs%eefld_wavelen
   write(*,9000) " ! -- Fock Bias Potentials: -- !"
   write(*,8180) inputs%fockbias_is_active, inputs%fockbias_is_shaped, &
                 inputs%fockbias_timeamp0
   write(*,8181) inputs%fockbias_timegrow, inputs%fockbias_timefall
   write(*,8182) inputs%fockbias_readfile

! General
9000 FORMAT(A)
8000 FORMAT(2x, "Natom = ", I6, ", Nsol = ", I8, ", charge = ", I5, &
            ", Nunp = ", I5, ", open = ", L2, ",")
8001 FORMAT(2x, "Nmax = ", I5, ", int_basis = ", L2, ", basis_set = ", A25, ",")
8002 FORMAT(2x, "fitting_set = ", A25, ", DIIS = ", L2, ", NDIIS = ", I3, ",")
8003 FORMAT(2x, "hybrid_converg = ", L2, ", Gold = ", F14.8, ", Told = ", &
            F14.8, ",")
8004 FORMAT(2x,"Etold = ", F14.8, ", good_cut = ", F14.8, ", Rmax = ", F14.8, &
            ",")
8005 FORMAT(2x,"RmaxS = ", F14.8, ", IExch = ", I5, ", IGrid = ", I3, &
            ", IGrid2 = ", I3, ",")
8006 FORMAT(2x, "PredCoef = ", L2, ", initial_guess = ", I3, ", DBug = ", L2, &
            ",")
8007 FORMAT(2x, "n_ghosts = ", I5)
! I/O Control
8020 FORMAT(2x, "verbose = ", I3, ", style = ", L2, ", timers = ", I3, &
            ", writeXYZ = ", L2, ", writeForces = ", L2, ",")
8021 FORMAT(2x, "dipole = ", L2, ", mulliken = ", L2, ", lowdin = ", L2, &
            ", fukui = ", L2, ", print_coeffs = ", L2, ",")
8022 FORMAT(2x, "VCInp = ", L2, ", FRestartIn = ", A25, ", restart_freq = ", &
            I5, ",")
8023 FORMAT(2x, "FRestart = ", A25, ", TDRestart = ", L2, ", writeDens = ", L2,&
            ",")
8024 FORMAT(2x, "TD_rst_freq = ", I6, ", gaussian_convert = ", L2,",")
8025 FORMAT(2x, "rst_dens = ", I3)
! TDDFT and Fields
8040 FORMAT(2x, "timeDep = ", I2, ", NTDStep = ", i10, ", TDStep = ", F14.8, &
           ", propagator = ", I2, ",")
8041 FORMAT(2x, "NBCH = ", I4, ", field = ", L2, ", a0 = ", F14.8, &
            ", epsilon = ", F14.8, ",")
8042 FORMAT(2x, "Fx = ", F14.8, ", Fy = ", F14.8, ", Fz = ", F14.8, &
            ", n_fields_iso = ", I5, ",")
8043 FORMAT(2x,"n_fields_aniso = ", I5, ", field_iso_file = ", A25, ",")
8044 FORMAT(2x,"field_aniso_file = ", A25, ", td_do_pop = ", I5)
! ECP
8060 FORMAT(2x, "ECPMode = ", L2, ", ECPTypes = ", I3, ", TipeECP = ", A25, ",")
8061 FORMAT(2x, "Fock_ECP_read = ", L2, ", Fock_ECP_write = ", L2, &
            ", cutECP = ", L2, ", cut2_0 = ", F14.8,",")
8062 FORMAT(2x, "cut3_0 = ", F14.8, ", verbose_ECP = ", I2, ", ECP_debug = ", &
            L2, ", fullTimer_ECP = ", L2, ",")
8063 FORMAT(2x, "local_nonlocal = ", I2, ", ECP_full_range_int = ", L2)
! Minimizations and restraints
8080 FORMAT(2x, "steep = ", L2, ", minimzation_steep = ", F14.8, &
            ", energy_cut = ", F14.8, ",")
8081 FORMAT(2x, "force_cut = ", F14.8, ", n_min_steeps = ", I5, &
            ", lineal_search = ", L2, ",")
8082 FORMAT(2x, "n_points = ", I5, ", number_restr = ", I5)
! CubeGen
8100 FORMAT(2x, "CubeGen_only = ", L2, ", cube_res = ", I5, ", cube_sel = ",   &
            I5, ", cube_dens = ", L2, ",")
8101 FORMAT(2x, "cube_orb = ", L2, ", cube_elec = ", L2, ", cube_dens_file = ",&
            A25, ",")
8102 FORMAT(2x, "cube_orb_file = ", A25, ", cube_elec_file = ", A25)
! GPU Options
8120 FORMAT(2x, "energy_all_iterations = ", L2, ", assign_all_functions = ", &
            L2, ", remove_zero_weights = ", L2, ",")
8121 FORMAT(2x, "max_function_exponent = ", I5, ", free_global_memory = ", &
            F14.8, ",")
8122 FORMAT(2x, "sphere_radius = ", F14.8, ", little_cube_size = ", F14.8, ",")
8123 FORMAT(2x, "min_points_per_cube = ", I5)
! Transport and DFTB
8140 FORMAT(2x, "transport_calc = ", L2, ", generate_rho0 = ", L2, &
            ", driving_rate = ", F14.8, ",")
8141 FORMAT(2x, "gate_field = ", L2, ", pop_drive = ", I3, &
            ", save_charge_freq = ", I5, ", DFTB_calc = ", L2, ",")
8142 FORMAT(2x, "MTB = ", I5, ", alfaTB = ", F14.8, ", betaTB = ", F14.8, ",")
8143 FORMAT(2x, "gammaTB = ", F14.8, ", VBias_TB = ", F14.8, ", start_TDTB = ",&
            I5, ",")
8144 FORMAT(2x, "end_TDTB = ", I5, ", end_BTB = ", I5, ", TBLoad = ", L2, &
            ", TBSave = ", L2)
8145 FORMAT(2x, "nbias=",I5)
! Ehrenfest
8160 FORMAT(2x, "ndyn_steps = ", I6, ", edyn_steps = ", I6, &
            ", nullify_forces = ", L2, ",")
8161 FORMAT(2x, "wdip_nfreq = ", I5, ", wdip_fname = ", A25, ", rsti_loads = ",&
            L2, ",")
8162 FORMAT(2x, "rsto_saves = ", L2, ", rsto_nfreq = ", L2, ", rsti_fname = ", &
            A25, ",")
8163 FORMAT(2x, "rsto_fname = ", A25, ", eefld_on = ", L2, ", eefld_ampx =", &
            F14.8, ",")
8164 FORMAT(2x, "eefld_ampy = ", F14.8, ", eefld_ampz = ", F14.8, ",")
8165 FORMAT(2x, "eefld_timeamp = ", F14.8, ", eefld_timepos = ", F14.8, ",")
8166 FORMAT(2x, "eefld_timegfh = ", L2, ", eefld_timegih = ", L2, &
            ", eefld_wavelen = ", F14.8)
! Fock Bias
8180 FORMAT(2x, "fockbias_is_active = ", L2, ", fockbias_is_shaped = ", L2, &
            ", fockbias_timeamp0 = ", F14.8, ",")
8181 FORMAT(2x, "fockbias_timegrow = ", F14.8, ", fockbias_timefall = ", F14.8,&
            ",")
8182 FORMAT(2x, "fockbias_readfile = ", A25)
   return
end subroutine lionml_write_dull


subroutine lionml_write_style()
   use lionml_data, only: lio_input_data, get_namelist
   type(lio_input_data) :: inputs

   call get_namelist(inputs)
   if (inputs%verbose .lt. 1) return

   ! LIO Header
   write(*,8000); write(*,8100); write(*,8001)

   ! General options and theory level
   write(*,8000); write(*,8101); write(*,8002)
   write(*,8200) inputs%natom         ; write(*,8201) inputs%nsol
   write(*,8202) inputs%charge        ; write(*,8203) inputs%Nunp
   write(*,8204) inputs%open          ; write(*,8205) inputs%nmax
   write(*,8206) inputs%int_basis     ; write(*,8207) inputs%basis_set
   write(*,8208) inputs%fitting_set   ; write(*,8209) inputs%diis
   write(*,8210) inputs%ndiis         ; write(*,8211) inputs%gold
   write(*,8212) inputs%told          ; write(*,8213) inputs%Etold
   write(*,8214) inputs%hybrid_converg; write(*,8215) inputs%good_cut
   write(*,8216) inputs%Rmax          ; write(*,8217) inputs%RmaxS
   write(*,8218) inputs%Iexch         ; write(*,8219) inputs%Igrid
   write(*,8220) inputs%Igrid2        ; write(*,8221) inputs%PredCoef
   write(*,8222) inputs%initial_guess ; write(*,8223) inputs%dbug
   write(*,8224) inputs%n_ghosts
   write(*,8003)

   ! File IO and Property calculations
   write(*,8000); write(*,8102); write(*,8002)
   write(*,8250) inputs%verbose     ; write(*,8251) inputs%style
   write(*,8252) inputs%timers      ; write(*,8253) inputs%writexyz
   write(*,8254) inputs%WriteForces ; write(*,8255) inputs%dipole
   write(*,8256) inputs%mulliken    ; write(*,8257) inputs%lowdin
   write(*,8258) inputs%fukui       ; write(*,8259) inputs%print_coeffs
   write(*,8260) inputs%restart_freq; write(*,8261) inputs%frestart
   write(*,8262) inputs%writedens   ; write(*,8263) inputs%td_rst_freq
   write(*,8264) inputs%vcinp       ; write(*,8265) inputs%Frestartin
   write(*,8266) inputs%Tdrestart   ; write(*,8267) inputs%gaussian_convert
   write(*,8268) inputs%rst_dens
   write(*,8003)

   ! TD-DFT and Fields
   write(*,8000); write(*,8103); write(*,8002)
   write(*,8300) inputs%timedep      ; write(*,8301) inputs%ntdstep
   write(*,8302) inputs%tdstep       ; write(*,8304) inputs%propagator
   write(*,8303) inputs%NBCH         ; write(*,8305) inputs%field
   write(*,8306) inputs%a0           ; write(*,8307) inputs%epsilon
   write(*,8308) inputs%Fx           ; write(*,8309) inputs%Fy
   write(*,8310) inputs%Fz           ; write(*,8311) inputs%nfields_iso
   write(*,8312) inputs%nfields_aniso; write(*,8313) inputs%field_iso_file
   write(*,8314) inputs%field_aniso_file
   write(*,8315) inputs%td_do_pop
   write(*,8003)

   ! Effective Core Potential
   write(*,8000); write(*,8104); write(*,8002)
   write(*,8350) inputs%Ecpmode       ; write(*,8351) inputs%Ecptypes
   write(*,8352) inputs%TipeECP
   call write_Zlist_ECP_style(inputs%ZlistECP, inputs%Ecptypes)
   write(*,8354) inputs%Fock_ECP_read ; write(*,8355) inputs%Fock_ECP_write
   write(*,8356) inputs%cutECP        ; write(*,8357) inputs%cut2_0
   write(*,8358) inputs%cut3_0        ; write(*,8359) inputs%Verbose_ECP
   write(*,8360) inputs%ECP_debug     ; write(*,8361) inputs%fulltimer_ECP
   write(*,8362) inputs%local_nonlocal; write(*,8363) inputs%ECP_full_range_int
   write(*,8003)

   ! Minimization and restraints
   write(*,8000); write(*,8105); write(*,8002)
   write(*,8370) inputs%steep       ; write(*,8371) inputs%minimzation_steep
   write(*,8372) inputs%Energy_cut  ; write(*,8373) inputs%Force_cut
   write(*,8374) inputs%n_min_steeps; write(*,8375) inputs%lineal_search
   write(*,8376) inputs%n_points    ; write(*,8377) inputs%number_restr
   write(*,8003)

   ! CUBEGEN
   write(*,8000); write(*,8106); write(*,8002)
   write(*,8400) inputs%Cubegen_only ; write(*,8401) inputs%Cube_Res
   write(*,8402) inputs%Cube_Dens    ; write(*,8403) inputs%Cube_Dens_file
   write(*,8404) inputs%Cube_Orb     ; write(*,8405) inputs%Cube_Sel
   write(*,8406) inputs%Cube_Orb_File; write(*,8407) inputs%Cube_Elec
   write(*,8408) inputs%Cube_Elec_File
   write(*,8003)

   ! GPU Options
   write(*,8000); write(*,8107); write(*,8002)
   write(*,8420) inputs%assign_all_functions
   write(*,8421) inputs%energy_all_iterations
   write(*,8422) inputs%remove_zero_weights
   write(*,8423) inputs%max_function_exponent
   write(*,8424) inputs%min_points_per_cube
   write(*,8425) inputs%little_cube_size
   write(*,8426) inputs%free_global_memory
   write(*,8427) inputs%sphere_radius
   write(*,8003)

   ! Transport and DFTB
   write(*,8000); write(*,8108); write(*,8002)
   write(*,8450) inputs%transport_calc; write(*,8451) inputs%generate_rho0
   write(*,8452) inputs%driving_rate  ; write(*,8453) inputs%gate_field
   write(*,8454) inputs%pop_drive     ; write(*,8455) inputs%save_charge_freq
   write(*,8456) inputs%dftb_calc     ; write(*,8457) inputs%MTB
   write(*,8458) inputs%alfaTB        ; write(*,8459) inputs%betaTB
   write(*,8460) inputs%gammaTB       ; write(*,8461) inputs%Vbias_TB
   write(*,8462) inputs%start_tdtb    ; write(*,8463) inputs%end_tdtb
   write(*,8464) inputs%end_bTB       ; write(*,8465) inputs%TBload
   write(*,8466) inputs%TBsave
   write(*,8003)

   ! Ehrenfest
   write(*,8000); write(*,8109); write(*,8002)
   write(*,8500) inputs%ndyn_steps    ; write(*,8501) inputs%edyn_steps
   write(*,8502) inputs%nullify_forces; write(*,8503) inputs%wdip_nfreq
   write(*,8504) inputs%wdip_fname    ; write(*,8505) inputs%rsti_loads
   write(*,8506) inputs%rsto_saves    ; write(*,8507) inputs%rsto_nfreq
   write(*,8508) inputs%rsti_fname    ; write(*,8509) inputs%rsto_fname
   write(*,8510) inputs%eefld_on      ; write(*,8511) inputs%eefld_ampx
   write(*,8512) inputs%eefld_ampy    ; write(*,8513) inputs%eefld_ampz
   write(*,8514) inputs%eefld_timeamp ; write(*,8515) inputs%eefld_timepos
   write(*,8516) inputs%eefld_timegfh ; write(*,8517) inputs%eefld_timegih
   write(*,8518) inputs%eefld_wavelen
   write(*,8003)

   ! Fock Bias Potentials
   write(*,8000); write(*,8110); write(*,8002)
   write(*,8550) inputs%fockbias_is_active
   write(*,8551) inputs%fockbias_is_shaped
   write(*,8552) inputs%fockbias_timeamp0
   write(*,8553) inputs%fockbias_timegrow
   write(*,8554) inputs%fockbias_timefall
   write(*,8555) inputs%fockbias_readfile
   write(*,8003)

   return;
8000 FORMAT(4x,"╔═════════════════════════════════", &
"═════════════════╗")
8001 FORMAT(4x,"╚═════════════════════════════════", &
"═════════════════╝")
8002 FORMAT(4x,"╠══════════════════════╦══════════", &
"═════════════════╣")
8003 FORMAT(4x,"╚══════════════════════╩══════════", &
"═════════════════╝")
8100 FORMAT(4x,"║                     LIO Input                    ║")
8101 FORMAT(4x,"║             General and Theory Level             ║")
8102 FORMAT(4x,"║             Input and Output Control             ║")
8103 FORMAT(4x,"║            RT-TDDFT and Field Options            ║")
8104 FORMAT(4x,"║             Effective Core Potential             ║")
8105 FORMAT(4x,"║            Minimization and Restraints           ║")
8106 FORMAT(4x,"║                     CubeGen                      ║")
8107 FORMAT(4x,"║                   GPU Options                    ║")
8108 FORMAT(4x,"║                Transport and DFTB                ║")
8109 FORMAT(4x,"║                Ehrenfest Dynamics                ║")
8110 FORMAT(4x,"║               Fock Bias Potentials               ║")

!System and Theory Level
8200 FORMAT(4x,"║  Natom               ║  ",17x,I6,2x,"║")
8201 FORMAT(4x,"║  Nsol                ║  ",15x,I8,2x,"║")
8202 FORMAT(4x,"║  Charge              ║  ",18x,I5,2x,"║")
8203 FORMAT(4x,"║  Nunp                ║  ",18x,I5,2x,"║")
8204 FORMAT(4x,"║  Open                ║  ",21x,L2,2x,"║")
8205 FORMAT(4x,"║  Nmax                ║  ",18x,I5,2x,"║")
8206 FORMAT(4x,"║  Int_Basis           ║  ",21x,L2,2x,"║")
8207 FORMAT(4x,"║  Basis_Set           ║  ",A25,"║")
8208 FORMAT(4x,"║  Fitting_Set         ║  ",A25,"║")
8209 FORMAT(4x,"║  Diis                ║  ",21x,L2,2x,"║")
8210 FORMAT(4x,"║  Ndiis               ║  ",20x,I3,2x,"║")
8211 FORMAT(4x,"║  Gold                ║  ",9x,F14.8,2x,"║")
8212 FORMAT(4x,"║  Told                ║  ",9x,F14.8,2x,"║")
8213 FORMAT(4x,"║  Etold               ║  ",9x,F14.8,2x,"║")
8214 FORMAT(4x,"║  Hybrid_converg      ║  ",21x,L2,2x,"║")
8215 FORMAT(4x,"║  Good_cut            ║  ",9x,F14.8,2x,"║")
8216 FORMAT(4x,"║  Rmax                ║  ",9x,F14.8,2x,"║")
8217 FORMAT(4x,"║  RmaxS               ║  ",9x,F14.8,2x,"║")
8218 FORMAT(4x,"║  Iexch               ║  ",18x,I5,2x,"║")
8219 FORMAT(4x,"║  Igrid               ║  ",20x,I3,2x,"║")
8220 FORMAT(4x,"║  Igrid2              ║  ",20x,I3,2x,"║")
8221 FORMAT(4x,"║  PredCoef            ║  ",21x,L2,2x,"║")
8222 FORMAT(4x,"║  Initial_guess       ║  ",18x,I5,2x,"║")
8223 FORMAT(4x,"║  Dbug                ║  ",21x,L2,2x,"║")
8224 FORMAT(4x,"║  n_ghosts            ║  ",18x,I5,2x,"║")
!IO Control
8250 FORMAT(4x,"║  Verbose             ║  ",20x,I3,2x,"║")
8251 FORMAT(4x,"║  Style               ║  ",21x,L2,2x,"║")
8252 FORMAT(4x,"║  Timers              ║  ",20x,I3,2x,"║")
8253 FORMAT(4x,"║  WriteXYZ            ║  ",21x,L2,2x,"║")
8254 FORMAT(4x,"║  WriteForces         ║  ",21x,L2,2x,"║")
8255 FORMAT(4x,"║  Dipole              ║  ",21x,L2,2x,"║")
8256 FORMAT(4x,"║  Mulliken            ║  ",21x,L2,2x,"║")
8257 FORMAT(4x,"║  Lowdin              ║  ",21x,L2,2x,"║")
8258 FORMAT(4x,"║  Fukui               ║  ",21x,L2,2x,"║")
8259 FORMAT(4x,"║  print_coeffs        ║  ",21x,L2,2x,"║")
8260 FORMAT(4x,"║  restart_freq        ║  ",17x,I6,2x,"║")
8261 FORMAT(4x,"║  Frestart            ║  ",A25,"║")
8262 FORMAT(4x,"║  WriteDens           ║  ",21x,L2,2x,"║")
8263 FORMAT(4x,"║  td_rst_freq         ║  ",17x,I6,2x,"║")
8264 FORMAT(4x,"║  VCinp               ║  ",21x,L2,2x,"║")
8265 FORMAT(4x,"║  Frestartin          ║  ",A25,"║")
8266 FORMAT(4x,"║  Tdrestart           ║  ",21x,L2,2x,"║")
8267 FORMAT(4x,"║  gaussian_convert    ║  ",21x,L2,2x,"║")
8268 FORMAT(4x,"║  rst_dens            ║  ",20x,I3,2x,"║")
! TD and Field options
8300 FORMAT(4x,"║  Timedep             ║  ",21x,I2,2x,"║")
8301 FORMAT(4x,"║  NTDstep             ║  ",13x,i10,2x,"║")
8302 FORMAT(4x,"║  TDstep              ║  ",9x,F14.8,2x,"║")
8303 FORMAT(4x,"║  Propagator          ║  ",21x,I2,2x,"║")
8304 FORMAT(4x,"║  NBCH                ║  ",19x,I4,2x,"║")
8305 FORMAT(4x,"║  Field               ║  ",21x,L2,2x,"║")
8306 FORMAT(4x,"║  A0                  ║  ",9x,F14.8,2x,"║")
8307 FORMAT(4x,"║  Epsilon             ║  ",9x,F14.8,2x,"║")
8308 FORMAT(4x,"║  Fx                  ║  ",9x,F14.8,2x,"║")
8309 FORMAT(4x,"║  Fy                  ║  ",9x,F14.8,2x,"║")
8310 FORMAT(4x,"║  Fz                  ║  ",9x,F14.8,2x,"║")
8311 FORMAT(4x,"║  n_fields_iso        ║  ",18x,I5,2x,"║")
8312 FORMAT(4x,"║  n_fields_aniso      ║  ",18x,I5,2x,"║")
8313 FORMAT(4x,"║  field_iso_file      ║  ",A25,"║")
8314 FORMAT(4x,"║  field_aniso_file    ║  ",A25,"║")
8315 FORMAT(4x,"║  td_do_pop           ║  ",18x,I5,2x,"║")
!Effective Core Potential
8350 FORMAT(4x,"║  Ecpmode             ║  ",21x,L2,2x,"║")
8351 FORMAT(4x,"║  Ecptypes            ║  ",20x,I3,2x,"║")
8352 FORMAT(4x,"║  TipeECP             ║  ",A25,"║")
8354 FORMAT(4x,"║  Fock_ECP_read       ║  ",21x,L2,2x,"║")
8355 FORMAT(4x,"║  Fock_ECP_write      ║  ",21x,L2,2x,"║")
8356 FORMAT(4x,"║  cutECP              ║  ",21x,L2,2x,"║")
8357 FORMAT(4x,"║  cut2_0              ║  ",9x,F14.8,2x,"║")
8358 FORMAT(4x,"║  cut3_0              ║  ",9x,F14.8,2x,"║")
8359 FORMAT(4x,"║  Verbose_ECP         ║  ",21x,I2,2x,"║")
8360 FORMAT(4x,"║  ECP_debug           ║  ",21x,L2,2x,"║")
8361 FORMAT(4x,"║  fulltimer_ECP       ║  ",21x,L2,2x,"║")
8362 FORMAT(4x,"║  local_nonlocal      ║  ",21x,I2,2x,"║")
8363 FORMAT(4x,"║  ECP_full_range_int  ║  ",21x,L2,2x,"║")
! Minimization and restraints
8370 FORMAT(4x,"║  Steep               ║  ",21x,L2,2x,"║")
8371 FORMAT(4x,"║  minimzation_steep   ║  ",9x,F14.8,2x,"║")
8372 FORMAT(4x,"║  Energy_cut          ║  ",9x,F14.8,2x,"║")
8373 FORMAT(4x,"║  Force_cut           ║  ",9x,F14.8,2x,"║")
8374 FORMAT(4x,"║  n_min_steeps        ║  ",18x,I5,2x,"║")
8375 FORMAT(4x,"║  lineal_search       ║  ",21x,L2,2x,"║")
8376 FORMAT(4x,"║  n_points            ║  ",18x,I5,2x,"║")
8377 FORMAT(4x,"║  number_restr        ║  ",18x,I5,2x,"║")
! Cubegen
8400 FORMAT(4x,"║  Cubegen_only        ║  ",21x,L2,2x,"║")
8401 FORMAT(4x,"║  Cube_Res            ║  ",18x,I5,2x,"║")
8402 FORMAT(4x,"║  Cube_Dens           ║  ",21x,L2,2x,"║")
8403 FORMAT(4x,"║  Cube_Dens_file      ║  ",A25,"║")
8404 FORMAT(4x,"║  Cube_Orb            ║  ",21x,L2,2x,"║")
8405 FORMAT(4x,"║  Cube_Sel            ║  ",18x,I5,2x,"║")
8406 FORMAT(4x,"║  Cube_Orb_File       ║  ",A25,"║")
8407 FORMAT(4x,"║  Cube_Elec           ║  ",21x,L2,2x,"║")
8408 FORMAT(4x,"║  Cube_Elec_File      ║  ",A25,"║")
! GPU options
8420 FORMAT(4x,"║  assign_all_functions║  ",21x,L2,2x,"║")
8421 FORMAT(4x,"║  energy_all_iteration║  ",21x,L2,2x,"║")
8422 FORMAT(4x,"║  remove_zero_weights ║  ",21x,L2,2x,"║")
8423 FORMAT(4x,"║  max_function_exponen║  ",18x,I5,2x,"║")
8424 FORMAT(4x,"║  min_points_per_cube ║  ",18x,I5,2x,"║")
8425 FORMAT(4x,"║  little_cube_size    ║  ",9x,F14.8,2x,"║")
8426 FORMAT(4x,"║  free_global_memory  ║  ",9x,F14.8,2x,"║")
8427 FORMAT(4x,"║  sphere_radius       ║  ",9x,F14.8,2x,"║")
! Transport and DFTB
8450 FORMAT(4x,"║  transport_calc      ║  ",21x,L2,2x,"║")
8451 FORMAT(4x,"║  generate_rho0       ║  ",21x,L2,2x,"║")
8452 FORMAT(4x,"║  driving_rate        ║  ",9x,F14.8,2x,"║")
8453 FORMAT(4x,"║  gate_field          ║  ",21x,L2,2x,"║")
8454 FORMAT(4x,"║  pop_drive           ║  ",20x,I3,2x,"║")
8455 FORMAT(4x,"║  save_charge_freq    ║  ",18x,I5,2x,"║")
8456 FORMAT(4x,"║  dftb_calc           ║  ",21x,L2,2x,"║")
8457 FORMAT(4x,"║  MTB                 ║  ",18x,I5,2x,"║")
8458 FORMAT(4x,"║  alfaTB              ║  ",9x,F14.8,2x,"║")
8459 FORMAT(4x,"║  betaTB              ║  ",9x,F14.8,2x,"║")
8460 FORMAT(4x,"║  gammaTB             ║  ",9x,F14.8,2x,"║")
8461 FORMAT(4x,"║  Vbias_TB            ║  ",9x,F14.8,2x,"║")
8462 FORMAT(4x,"║  start_tdtb          ║  ",18x,I5,2x,"║")
8463 FORMAT(4x,"║  end_tdtb            ║  ",18x,I5,2x,"║")
8464 FORMAT(4x,"║  end_bTB             ║  ",18x,I5,2x,"║")
8465 FORMAT(4x,"║  TBload              ║  ",21x,L2,2x,"║")
8466 FORMAT(4x,"║  TBsave              ║  ",21x,L2,2x,"║")
! Ehrenfest
8500 FORMAT(4x,"║  ndyn_steps          ║  ",17x,I6,2x,"║")
8501 FORMAT(4x,"║  edyn_steps          ║  ",17x,I6,2x,"║")
8502 FORMAT(4x,"║  nullify_forces      ║  ",21x,L2,2x,"║")
8503 FORMAT(4x,"║  wdip_nfreq          ║  ",18x,I5,2x,"║")
8504 FORMAT(4x,"║  wdip_fname          ║  ",A25,"║")
8505 FORMAT(4x,"║  rsti_loads          ║  ",21x,L2,2x,"║")
8506 FORMAT(4x,"║  rsto_saves          ║  ",21x,L2,2x,"║")
8507 FORMAT(4x,"║  rsto_nfreq          ║  ",21x,L2,2x,"║")
8508 FORMAT(4x,"║  rsti_fname          ║  ",A25,"║")
8509 FORMAT(4x,"║  rsto_fname          ║  ",A25,"║")
8510 FORMAT(4x,"║  eefld_on            ║  ",21x,L2,2x,"║")
8511 FORMAT(4x,"║  eefld_ampx          ║  ",9x,F14.8,2x,"║")
8512 FORMAT(4x,"║  eefld_ampy          ║  ",9x,F14.8,2x,"║")
8513 FORMAT(4x,"║  eefld_ampz          ║  ",9x,F14.8,2x,"║")
8514 FORMAT(4x,"║  eefld_timeamp       ║  ",9x,F14.8,2x,"║")
8515 FORMAT(4x,"║  eefld_timepos       ║  ",9x,F14.8,2x,"║")
8516 FORMAT(4x,"║  eefld_timegfh       ║  ",21x,L2,2x,"║")
8517 FORMAT(4x,"║  eefld_timegih       ║  ",21x,L2,2x,"║")
8518 FORMAT(4x,"║  eefld_wavelen       ║  ",9x,F14.8,2x,"║")
! Fock Bias Potentials
8550 FORMAT(4x,"║  fockbias_is_active  ║  ",21x,L2,2x,"║")
8551 FORMAT(4x,"║  fockbias_is_shaped  ║  ",21x,L2,2x,"║")
8552 FORMAT(4x,"║  fockbias_timeamp0   ║  ",9x,F14.8,2x,"║")
8553 FORMAT(4x,"║  fockbias_timegrow   ║  ",9x,F14.8,2x,"║")
8554 FORMAT(4x,"║  fockbias_timefall   ║  ",9x,F14.8,2x,"║")
8555 FORMAT(4x,"║  fockbias_readfile   ║  ",A25,"║")
end subroutine lionml_write_style

subroutine write_Zlist_ECP_dull(ZlistECP, D)
   implicit none
   integer, intent(in) :: ZlistECP(128)
   integer, intent(in) :: D
   integer :: icount, kcount, lines, index

   write(*, 9000, advance='no') "  ZListECP = "
   if (D .lt. 1) then
      write(*,*)
      return;
   else
      if (D .lt. 21) then
         do icount = 1, D-1
            write(*, 9002, advance='no') ZlistECP(icount)
         enddo
      else
         do icount = 1, 19
            write(*, 9002, advance='no') ZlistECP(icount)
         enddo
         write(*, 9002) ZlistECP(20)
         write(*, 9001, advance='no') "  "
         if (D .lt. 46) then
            do icount = 21, D-1
               write(*, 9002, advance='no') ZlistECP(icount)
            enddo
         else
            write(*, 9001, advance='no') "  "
            lines = (D - 20) / 25
            do icount=1, lines - 1
               do kcount=1, 24
                  index = 20 + kcount*icount
                  write(*, 9002, advance='no') ZlistECP(index)
               enddo
               index = 20 + 25*icount
               write(*, 9002) ZlistECP(index)
               write(*, 9001, advance='no') "  "
            enddo
            do kcount = D - (20 + 25*(lines - 1)), D-1
               index = 20 + kcount*lines
               write(*, 9002, advance='no') ZlistECP(index)
            enddo
         endif
      endif
   endif
   write(*,9002) ZlistECP(D)

   return;

9000 FORMAT(A)
9001 FORMAT(A2)
9002 FORMAT(I3)
end subroutine write_Zlist_ECP_dull

subroutine write_Zlist_ECP_style(ZlistECP, D)
   implicit none
   integer, intent(in) :: ZlistECP(128)
   integer, intent(in) :: D
   integer :: icount, kcount, lines, rest

   if (D .lt. 6) then
     if (D .eq. 1) write(*,9538) ZlistECP(1)
     if (D .eq. 2) write(*,9539) ZlistECP(1:2)
     if (D .eq. 3) write(*,9540) ZlistECP(1:3)
     if (D .eq. 4) write(*,9541) ZlistECP(1:4)
     if (D .eq. 5) write(*,9542) ZlistECP(1:5)
   else
      lines = D / 6
      rest  = mod(D, 6)
      write(*,9543) ZlistECP(1:6)
      do icount = 1, lines-1
         kcount = 6*icount + 1
         write(*,9544) ZlistECP(kcount:kcount+5)
      enddo
      if (rest .eq. 1) write(*,9545) ZlistECP(6*lines+1:D)
      if (rest .eq. 2) write(*,9546) ZlistECP(6*lines+1:D)
      if (rest .eq. 3) write(*,9547) ZlistECP(6*lines+1:D)
      if (rest .eq. 4) write(*,9548) ZlistECP(6*lines+1:D)
      if (rest .eq. 5) write(*,9549) ZlistECP(6*lines+1:D)
   endif

   return;
9538 FORMAT(4x,"║  Zlistecp            ║ ",I3,"                       ║")
9539 FORMAT(4x,"║  Zlistecp            ║ ",I3,I3,"                    ║")
9540 FORMAT(4x,"║  Zlistecp            ║ ",I3,I3,I3,"                 ║")
9541 FORMAT(4x,"║  Zlistecp            ║ ",I3,I3,I3,I3,"              ║")
9542 FORMAT(4x,"║  Zlistecp            ║ ",I3,I3,I3,I3,I3,"           ║")
9543 FORMAT(4x,"║  Zlistecp            ║ ",I3,I3,I3,I3,I3,I3,"        ║")
9544 FORMAT(4x,"║                      ║ ",I3,I3,I3,I3,I3,I3,"        ║")
9545 FORMAT(4x,"║                      ║ ",I3"                        ║")
9546 FORMAT(4x,"║                      ║ ",I3,I3,"                    ║")
9547 FORMAT(4x,"║                      ║ ",I3,I3,I3,"                 ║")
9548 FORMAT(4x,"║                      ║ ",I3,I3,I3,I3"               ║")
9549 FORMAT(4x,"║                      ║ ",I3,I3,I3,I3,I3"            ║")
end subroutine write_Zlist_ECP_style
end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
