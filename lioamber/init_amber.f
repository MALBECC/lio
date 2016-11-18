      subroutine init_lio_amber(natomin,Izin,nclatom,charge
     > , basis_i, output_i, fcoord_i, fmulliken_i, frestart_i
     > , frestartin_i, verbose_i, OPEN_i, NMAX_i, NUNP_i, VCINP_i
     > , GOLD_i, told_i, rmax_i, rmaxs_i, predcoef_i, idip_i
     > , writexyz_i, intsoldouble_i, DIIS_i, ndiis_i, dgtrig_i, Iexch_i
     > , integ_i, DENS_i , IGRID_i, IGRID2_i , timedep_i , tdstep_i 
     > , ntdstep_i, field_i, exter_i, a0_i, epsilon_i, Fx_i
     > , Fy_i, Fz_i, NBCH_i, propagator_i, writedens_i, tdrestart_i
#ifdef MOD_AMBER
     > , basis_set_i, fitting_set_i, int_basis_i
     > , cubegen_only_i, cuberes_i
     > , cubedens_i, cubedensfile_i
     > , cubeorb_i, cubesel_i, cubeorbfile_i
     > , restart_freq_i, energy_freq_i)
#else
     > )
#endif

      use garcha_mod
      use ECP_mod, only : ecpmode, ecptypes, tipeECP, ZlistECP
     > ,cutECP,local_nonlocal, ecp_debug,ecp_full_range_int
     > ,verbose_ECP,Cnorm,FOCK_ECP_read, FOCK_ECP_write,Fulltimer_ECP
     > ,cut2_0,cut3_0

c      use qmmm_module, only : qmmm_struct,qmmm_nml
      implicit real*8 (a-h,o-z)
c
c PARAMETERS - DYNAMICAL VECTOR ONLY -------------------
c
c ngDyn : number of atoms * number basis functions
c ngdDyn: number of atoms * number auxiliar functions
c Ngrid : number of grid points (LS-SCF part)
c norbit : number of MO
c
c Ngrid may be set to 0 , in the case of using Num. Integ.
c
c      parameter (ngDyn=700)
c      parameter (ngdDyn=850)


c      include 'param'
        integer , intent(in) :: charge, nclatom
       integer , intent(in)  :: natomin
         integer , intent(in)  :: Izin(natomin)
       character(len=20) :: basis_i
#ifdef MOD_AMBER
       character(len=40) :: basis_set_i
       character(len=40) :: fitting_set_i
       logical :: int_basis_i,cubegen_only_i,cubedens_i,cubeorb_i
       integer :: cuberes_i, cubesel_i
       character(len=20) :: cubedensfile_i,cubeorbfile_i
       integer :: restart_freq_i, energy_freq_i
#endif
       character(len=20) :: output_i
       character(len=20) :: fcoord_i
         character(len=20) :: fmulliken_i
       character(len=20) :: frestart_i
       character(len=20) :: frestartin_i
       logical :: verbose_i
       logical :: OPEN_i
       integer :: NMAX_i
       integer :: NUNP_i
       logical :: VCINP_i
         real*8  :: GOLD_i
       real*8  :: told_i
       real*8  :: rmax_i
       real*8  :: rmaxs_i
       logical :: predcoef_i
       integer :: idip_i
       logical :: writexyz_i
       logical :: intsoldouble_i
       logical :: DIIS_i
       integer :: ndiis_i
       real*8  :: dgtrig_i
       integer :: Iexch_i
       logical :: integ_i
       logical :: DENS_i
       integer :: IGRID_i
       integer :: IGRID2_i
       integer :: timedep_i
       logical :: field_i
       logical :: exter_i
       real*8  :: tdstep_i
       integer  :: ntdstep_i
       real*8  :: a0_i
       real*8  :: epsilon_i
       real*8  :: Fx_i
       real*8  :: Fy_i
       real*8  :: Fz_i
       integer  :: NBCH_i
       integer :: propagator_i
       logical :: writedens_i
       logical :: tdrestart_i
       
       call lio_defaults()  

       basis= basis_i
#ifdef MOD_AMBER
       basis_set=basis_set_i
       fitting_set=fitting_set_i
       int_basis=int_basis_i
       cubegen_only = cubegen_only_i
       cube_res = cuberes_i
       cube_dens = cubedens_i
       cube_dens_file = cubedensfile_i
       cube_orb = cubeorb_i
       cube_sel = cubesel_i
       cube_orb_file = cubeorbfile_i
       restart_freq = restart_freq_i
       energy_freq = energy_freq_i
#else
       basis_set="DZVP"
       fitting_set="DZVP Coulomb Fitting"
       int_basis=.false.
       cubegen_only = .false.
       cube_res = 40
       cube_dens =.false.
       cube_dens_file = 'dens.cube'
       cube_orb = .false.
       cube_sel = 0
       cube_orb_file = "orb.cube"
       restart_freq = 1
       energy_freq = 1
#endif
       Output= output_i
       fcoord=fcoord_i
       fmulliken=fmulliken_i
       frestart= frestart_i
       frestartin=frestartin_i
        verbose = verbose_i
       OPEN=OPEN_i
       NMAX= NMAX_i
       NUNP=NUNP_i
       VCINP=VCINP_i
       GOLD=GOLD_i
       told=told_i
       rmax=rmax_i
       rmaxs=rmaxs_i
       predcoef=predcoef_i
       idip=idip_i
       writexyz=writexyz_i
       intsoldouble= intsoldouble_i
       DIIS=DIIS_i
       ndiis=ndiis_i
       dgtrig= dgtrig_i
       Iexch=Iexch_i
       integ=integ_i
       DENS=DENS_i
       IGRID=IGRID_i
       IGRID2=IGRID2_i
       timedep=timedep_i
       field=field_i
       exter=exter_i
       tdstep=tdstep_i
       ntdstep= ntdstep_i
       a0=a0_i
       epsilon=epsilon_i
       Fx=Fx_i
       Fy=Fy_i
       Fz=Fz_i
       NBCH=NBCH_i
       propagator=propagator_i
       writedens=writedens_i
       tdrestart=tdrestart_i

      call lio_init() 
      end
