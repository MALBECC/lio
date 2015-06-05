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


c      parameter (norbit=800,Ngrid=0)


         natom=natomin

c       integer, intent(in) :: Iiiz(natom)
         ntatom=natom+nclatom
         ntatom=ntatom ! the number of clasical atoms can change
         ngnu=natom*ng0
         ngdnu=natom*ngd0
         ngDyn=ngnu
         ngdDyn=ngdnu
c
        ng3=4*ngDyn
c para version en memoria
      ng2=5*ngDyn*(ngDyn+1)/2+3*ngdDyn*(ngdDyn+1)/2+
     >           ngDyn+ngDyn*norbit+Ngrid

c      write(*,*) 'ng2 en init',ng2,ngDyn,ngdDyn,norbit,Ngrid

      allocate(X(ngDyn,ng3),XX(ngdDyn,ngdDyn))
      allocate(RMM(ng2),RMM1(ng2),RMM2(ng2), RMM3(ng2))
      
       allocate (c(ngnu,nl),a(ngnu,nl),Nuc(ngnu),ncont(ngnu)
     >  ,cx(ngdnu,nl),ax(ngdnu,nl),Nucx(ngdnu),ncontx(ngdnu)
     > ,cd(ngdnu,nl),ad(ngdnu,nl),Nucd(ngdnu),ncontd(ngdnu)
     > ,indexii(ngnu),indexiid(ngdnu))

      allocate (r(ntatom,3),v(ntatom,3),rqm(natom,3),Em(ntatom)
     >,Rm(ntatom),pc(ntatom),Iz(natom),nnat(ntatom),af(natom*ngd0),
     >  B(natom*ngd0,3))
      allocate(d(natom,natom))
         Iz=Izin
      if(verbose) then
      write(6,*) '---------Lio options-------'
      write(6,*)      '  OPEN ', OPEN
       write(6,*)     '  NMAX ', NMAX
       write(6,*)     '  NUNP ', NUNP
       write(6,*)     '  VCINP ', VCINP
!   write(6,*)     '  GOLD ', GOLD
       write(6,*)     '  told ', told
!   write(6,*)     '  rmax ', rmax
!   write(6,*)     '  rmaxs ', rmaxs
!     write(6,*)     '  predcoef ', predcoef
!     write(6,*)     '  idip ', idip
       write(6,*)     '  writexyz ', writexyz
       write(6,*)     '  DIIS ', DIIS
       write(6,*)     '  ndiis ', ndiis
       write(6,*)     '  Iexch ', Iexch
!   write(6,*)     '  integ ', integ
!     write(6,*)     '  DENS ' ,  DENS
       write(6,*)     '  IGRID ', IGRID
       write(6,*)     '  IGRID2 ', IGRID2
       write(6,*)     '  timedep ', timedep
       write(6,*)     '  tdstep ', tdstep
       write(6,*)     '  ntdstep ', ntdstep
       write(6,*)     '  field ', field
       write(6,*)     '  exter ', exter
       write(6,*)     '  a0 ', a0
       write(6,*)     '  epsilon ', epsilon
       write(6,*)     '  Fx ', Fx
       write(6,*)     '  Fy ', Fy
       write(6,*)     '  Fz ', Fz
       write(6,*)     '  NBCH ', NBCH
       write(6,*)     '  propagator ', propagator
       write(6,*)     '  writedens ', writedens
       write(6,*)     '  tdrestart ', tdrestart
        write(6,*) '-----end Lio options-------'
       endif


c        write(92,*) izin
c      write(*,*) (ngDyn*ng3+ngdDyn**2+ng2+ngDyn*
c     > (NgDyn+1)/2*NgdDyn)*8*1.0D-06, '  Memoria en MB'
c      write(*,*)   ng2
!#ifdef G2G
      call g2g_init()
!#endif
        nqnuc=0
       do i=1,natom
        nqnuc=nqnuc+Iz(i)
        enddo

c !!! REVISAR QUE nco se un numero ENTERO ¡¡¡¡¡¡

        nco=((nqnuc - charge)-Nunp)/2

c
c        write(*,*) 'NCO=',NCO
c       write(*,*) natom,ntatom,ngDyn,ngdDyn,ng0,ngd0
c       write(*,*) ng2,ngDyn,ngdDyn
c--------------------------------------------------------
      call drive(ng2,ngDyn,ngdDyn)
c        write(*,*) 'Lio init amber'

      end
c---------------------------------------------------------------------
