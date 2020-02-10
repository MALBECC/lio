#include "../include/dprec.fh"
module qm2_extern_lio_module
! ----------------------------------------------------------------
! Interface for lio based QM MD 
!
! Currently supports:
! pure QM and QM/MM
!
! Initial implementation by
! Matthew Clark
! under supervision of
! Andreas Goetz and Ross Walker (SDSC) 
! 
! Date: February 2011
!
!
! Modified by Mariano C. González Lebrero
! Date: March 2013 march
!
! ----------------------------------------------------------------
  implicit none

  private
  public :: get_lio_forces

  type lio_nml_type
     character(len=20) :: basis
     character(len=20) :: output
     character(len=20) :: fcoord
     character(len=20) :: fmulliken
     character(len=20) :: frestart
     character(len=20) :: frestartin
     logical :: verbose
     logical :: OPEN
     integer :: NMAX
     integer :: NUNP
     logical :: VCINP
     real*8  :: GOLD
     real*8  :: told
     real*8  :: rmax
     real*8  :: rmaxs
     logical :: predcoef
     integer :: idip
     logical :: writexyz
     logical :: intsoldouble
     logical :: DIIS
     integer :: ndiis
     real*8  :: dgtrig
     integer :: Iexch
     logical :: integ
     logical :: DENS
     integer :: IGRID 
     integer :: IGRID2 
     integer :: timedep 
     real*8  :: tdstep 
     integer  :: ntdstep 
     logical :: field
     logical :: exter 
     real*8  :: a0
     real*8  :: epsilon
     real*8  :: Fx
     real*8  :: Fy
     real*8  :: Fz
     integer  :: NBCH
     integer :: propagator
     logical :: writedens
     logical :: tdrestart

  end type lio_nml_type

contains

  ! --------------------------------------------
  ! Get QM energy and forces from Lio 
  ! --------------------------------------------
  subroutine get_lio_forces( do_grad, nstep, ntpr_default, id, nqmatoms, qmcoords,&
    nclatoms, clcoords, escf, dxyzqm, dxyzcl)

    use memory_module, only:x,lvel
    use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS, ZERO
    use file_io_dat
    use qmmm_module, only : qmmm_struct, qmmm_nml
!    use garcha_mod


    implicit none

#  include "../include/md.h" ! FFR-Ehrenfest


    logical, intent(in) :: do_grad              ! Return gradient/not
    integer, intent(in) :: nstep                ! MD step number
    integer, intent(in) :: ntpr_default         ! frequency of printing
    character(len=3), intent(in) :: id          ! ID number for PIMD or REMD
    integer, intent(in) :: nqmatoms             ! Number of QM atoms
    _REAL_,  intent(in) :: qmcoords(3,nqmatoms) ! QM atom coordinates
    _REAL_              :: qmvels(3,nqmatoms)   ! QM atom velocities (of previous step) (FFR-Ehrenfest)
    integer, intent(in) :: nclatoms             ! Number of MM atoms
    _REAL_,  intent(in) :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
    _REAL_, intent(out) :: escf                 ! SCF energy
    _REAL_, intent(out) :: dxyzqm(3,nqmatoms)   ! SCF QM force
    _REAL_, intent(out) :: dxyzcl(3,nclatoms)   ! SCF MM force
    _REAL_              :: dipxyz(3), dipole    ! Dipole moment

    type(lio_nml_type), save     :: lio_nml
    logical, save                :: first_call = .true.
    integer                      :: nn,i        ! FFR-Ehrenfest: added nn
    integer                      :: printed =-1 ! Used to tell if we have printed this step yet 
                                                ! since the same step may be called multiple times
    ! for system call
    integer :: system
    integer :: stat


    ! Setup on first call
    if ( first_call ) then
       first_call = .false.
       write (6,'(/,a,/)') '  >>> Running calculations with Lio <<<'
       call get_namelist_lio(lio_nml)
       call print_namelist(lio_nml)
!------------------------------------------------------------------------------!
! FFR - EHRENFEST ADAPTATION
!------------------------------------------------------------------------------!
       write (6,'(/,a,/)') '  >>> Using code in /LioDevelop/Amber14eh2 <<<'
!       call init_lio_amber(nqmatoms,qmmm_struct%iqm_atomic_numbers,nclatoms, &
       call init_lioamber_ehren(nqmatoms,qmmm_struct%iqm_atomic_numbers,nclatoms, &
          qmmm_nml%qmcharge, lio_nml%basis, lio_nml%output, lio_nml%fcoord, &
          lio_nml%fmulliken, lio_nml%frestart, lio_nml%frestartin, &
          lio_nml%verbose, lio_nml%OPEN, lio_nml%NMAX, lio_nml%NUNP, &
          lio_nml%VCINP, lio_nml%GOLD, lio_nml%told, lio_nml%rmax, &
          lio_nml%rmaxs, lio_nml%predcoef, lio_nml%idip, lio_nml%writexyz, &
          lio_nml%intsoldouble, lio_nml%DIIS, lio_nml%ndiis, lio_nml%dgtrig, &
          lio_nml%Iexch, lio_nml%integ, lio_nml%DENS, lio_nml%IGRID, &
          lio_nml%IGRID2, lio_nml%timedep, lio_nml%tdstep, lio_nml%ntdstep, &
          lio_nml%field, lio_nml%exter, lio_nml%a0, lio_nml%epsilon, &
          lio_nml%Fx, lio_nml%Fy, lio_nml%Fz, lio_nml%NBCH, &
!          lio_nml%propagator, lio_nml%writedens, lio_nml%tdrestart)
          lio_nml%propagator, lio_nml%writedens, lio_nml%tdrestart,dt)
!------------------------------------------------------------------------------!
    end if

!           if(writexyz) then
!            do i=1,nqmatoms
!            jj=(qmmm_nml%iqmatoms(i)-1)*3
!       write(18,345) qmmm_struct%iqm_atomic_numbers(i),coords(jj+1),&
!           coords(jj+2),coords(jj+3)
!         enddo
!!           endif


!!      write(665,123) x(lvel:lvel+3*nqmatoms-1)*20.455d0
!!123 FORMAT(3(3X,f10.6))
!------------------------------------------------------------------------------!
! FFR - EHRENFEST ADAPTATION
!------------------------------------------------------------------------------!
      do nn=1,nqmatoms
        qmvels(1,nn)=x(lvel+3*qmmm_struct%iqmatoms(nn)-3)
        qmvels(2,nn)=x(lvel+3*qmmm_struct%iqmatoms(nn)-2)
        qmvels(3,nn)=x(lvel+3*qmmm_struct%iqmatoms(nn)-1)
      enddo

!      call SCF_in(escf,qmcoords,clcoords,nclatoms,dipxyz)
      call ehren_in( qmcoords, qmvels, clcoords, nclatoms, dipxyz, escf)
!------------------------------------------------------------------------------!


          escf=escf*CODATA08_AU_TO_KCAL

!      if (do_grad) then
         dxyzqm=0
         dxyzcl=0

      call dft_get_qm_forces(dxyzqm)
!      open(unit=44,file='dxyzqm')
!      write(44,*) dxyzqm  
      call dft_get_mm_forces(dxyzcl,dxyzqm)

       dxyzqm=dxyzqm*CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS


       dxyzcl=dxyzcl*CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS

!      open(unit=45,file='dxyzqm2')
!      open(unit=46,file='dxyzcl')
!      write(45,*) dxyzqm  
!      write(46,*) dxyzcl  

!      endif


!    escf = escf * CODATA08_AU_TO_KCAL

  end subroutine get_lio_forces

  ! ---------------------------------------------
  ! Read Lio namelist values from file mdin,
  ! use default values if none are present.
  ! ---------------------------------------------
    
 
  subroutine get_namelist_lio(lio_nml)
!    use garcha_mod
    implicit none
    type(lio_nml_type), intent(out) :: lio_nml

     character(len=20) :: basis
     character(len=20) :: output
     character(len=20) :: fcoord
     character(len=20) :: fmulliken
     character(len=20) :: frestart
     character(len=20) :: frestartin
     logical :: verbose
     logical :: OPEN
     integer :: NMAX
     integer :: NUNP
     logical :: VCINP
     real*8  :: GOLD
     real*8  :: told
     real*8  :: rmax
     real*8  :: rmaxs
     logical :: predcoef
     integer :: idip
     logical :: writexyz
     logical :: intsoldouble
     logical :: DIIS
     integer :: ndiis
     real*8  :: dgtrig
     integer :: Iexch
     logical :: integ
     logical :: DENS
     integer :: IGRID
     integer :: IGRID2
     integer :: timedep
     logical :: field
     logical :: exter
     real*8  :: tdstep
     integer  :: ntdstep
     real*8  :: a0
     real*8  :: epsilon
     real*8  :: Fx
     real*8  :: Fy
     real*8  :: Fz
     integer  :: NBCH
     integer :: propagator
     logical :: writedens
     logical :: tdrestart

     namelist /lio/ basis,output,fcoord,fmulliken,frestart,frestartin &
     ,verbose,OPEN,NMAX,Nunp,VCINP, &
     GOLD,told,rmax,rmaxs,predcoef, &
     idip,writexyz,intsoldouble,DIIS,ndiis,dgtrig, &
     Iexch,integ,dens,igrid,igrid2,timedep, tdstep, ntdstep, &
     propagator,NBCH, writedens, tdrestart, &
     field,exter,a0,epsilon,Fx,Fy,Fz

    integer :: ifind, ierr
    !defaults

     basis='basis'  ! name of the base file
     output='output'
     fcoord='qm.xyz'
     fmulliken='mulliken'
     frestart='restart.out'
     frestartin='restart.in' 
     verbose=.false.
     OPEN=.false.
     NMAX= 100
     NUNP= 0
     VCINP= .false.
     GOLD= 10.
     told=1.0D-6
     rmax=16
     rmaxs=5
     predcoef=.false.
     idip=1 
     writexyz=.true.
     intsoldouble=.true.
     DIIS=.true.
     ndiis=30
     dgtrig=100.
     Iexch=9
     integ=.true.
     DENS = .true.
     IGRID = 2
     IGRID2 = 2
     timedep = 0
     tdstep = 2.D-3
     ntdstep= 1
     field=.false.
     exter=.false.
     a0=1000.0
     epsilon=1.D0
     Fx=0.05
     Fy=0.05
     Fz=0.05
     NBCH=10
     propagator=1
     writedens=.true.
     tdrestart=.false.
    ! Read namelist
    rewind 5
    read(5,nml=lio,iostat=ierr)

    if ( ierr > 0 ) then
       call sander_bomb('get_namelist_lio (qm2_extern_lio_module)', &
            '&lio namelist read error', &
            'Please check your input.')
    else if ( ierr < 0 ) then
       write(6,'(a)') '&lio namelist read success'
    end if



     lio_nml%basis=basis
     lio_nml%output=output
     lio_nml%fcoord=fcoord
     lio_nml%fmulliken=fmulliken
     lio_nml%frestart=frestart
     lio_nml%frestartin=frestartin
     lio_nml%verbose=verbose
     lio_nml%OPEN=OPEN
     lio_nml%NMAX=NMAX
     lio_nml%NUNP=NUNP
     lio_nml%VCINP=VCINP
     lio_nml%GOLD=GOLD
     lio_nml%told=told
     lio_nml%rmax=rmax
     lio_nml%rmaxs=rmaxs
     lio_nml%predcoef=predcoef
     lio_nml%idip=idip
     lio_nml%writexyz=writexyz
     lio_nml%intsoldouble=intsoldouble
     lio_nml%DIIS=DIIS
     lio_nml%ndiis=ndiis
     lio_nml%dgtrig=dgtrig
     lio_nml%Iexch=Iexch
     lio_nml%integ=integ
     lio_nml%DENS =DENS 
     lio_nml%IGRID =IGRID 
     lio_nml%IGRID2 =IGRID2 
     lio_nml%timedep =timedep 
     lio_nml%tdstep = tdstep 
     lio_nml%ntdstep = ntdstep 
     lio_nml%field=field
     lio_nml%exter=exter
     lio_nml%a0=a0
     lio_nml%epsilon=epsilon
     lio_nml%Fx=Fx
     lio_nml%Fy=Fy
     lio_nml%Fz=Fz
     lio_nml%NBCH=NBCH
     lio_nml%propagator=propagator
     lio_nml%writedens=writedens
     lio_nml%tdrestart=tdrestart

  end subroutine get_namelist_lio

    subroutine print_namelist(lio_nml)
    implicit none
    type(lio_nml_type), intent(in) :: lio_nml
!     use garcha_mod
  write(6,*) '---------Lio options-------'
  write(6,*)      '  OPEN ', lio_nml%OPEN
   write(6,*)     '  NMAX ', lio_nml%NMAX
   write(6,*)     '  NUNP ', lio_nml%NUNP
   write(6,*)     '  VCINP ', lio_nml%VCINP
!   write(6,*)     '  GOLD ', GOLD
   write(6,*)     '  told ', lio_nml%told
!   write(6,*)     '  rmax ', rmax
!   write(6,*)     '  rmaxs ', rmaxs
!   write(6,*)     '  predcoef ', predcoef
!   write(6,*)     '  idip ', idip
   write(6,*)     '  writexyz ', lio_nml%writexyz
   write(6,*)     '  DIIS ', lio_nml%DIIS
   write(6,*)     '  ndiis ', lio_nml%ndiis
   write(6,*)     '  Iexch ', lio_nml%Iexch
!   write(6,*)     '  integ ', integ
!   write(6,*)     '  DENS ' ,  DENS
   write(6,*)     '  IGRID ', lio_nml%IGRID
   write(6,*)     '  IGRID2 ', lio_nml%IGRID2
   write(6,*)     '  timedep ', lio_nml%timedep
    if(lio_nml%timedep.gt.0) then
   write(6,*)     'Runing TDDFT calculation with Lio'
   write(6,*)     '  tdstep ', lio_nml%tdstep
   write(6,*)     '  ntdstep ', lio_nml%ntdstep
   write(6,*)     '  field ', lio_nml%field
   write(6,*)     '  exter ', lio_nml%exter
   write(6,*)     '  a0 ', lio_nml%a0
   write(6,*)     '  epsilon ', lio_nml%epsilon
   write(6,*)     '  Fx ', lio_nml%Fx 
   write(6,*)     '  Fy ', lio_nml%Fy
   write(6,*)     '  Fz ', lio_nml%Fz
   write(6,*)     '  NBCH ', lio_nml%NBCH 
   write(6,*)     '  propagator ', lio_nml%propagator
   write(6,*)     '  writedens ', lio_nml%writedens
   write(6,*)     '  writedens ', lio_nml%tdrestart
   endif
  write(6,*) '-----end Lio options-------'

     end subroutine print_namelist
!

end module qm2_extern_lio_module
