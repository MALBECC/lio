!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! MODULES
!
! Add line 4 starting with "use memory_module"
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine get_lio_forces( do_grad, nstep, ntpr_default, id, nqmatoms, qmcoords,&
       nclatoms, clcoords, escf, dxyzqm, dxyzcl)
  
       use memory_module, only:x,lvel
       use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS, ZERO
       use file_io_dat
       use qmmm_module, only : qmmm_struct, qmmm_nml
!      use garcha_mod
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! DECLARATIONS
!
! Add lines 3 (include...), 11 (qmvels) and 21 (nn declaration).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       implicit none

#      include "../include/md.h"

       logical, intent(in) :: do_grad              ! Return gradient/not
       integer, intent(in) :: nstep                ! MD step number
       integer, intent(in) :: ntpr_default         ! frequency of printing
       character(len=3), intent(in) :: id          ! ID number for PIMD or REMD
       integer, intent(in) :: nqmatoms             ! Number of QM atoms
       _REAL_,  intent(in) :: qmcoords(3,nqmatoms) ! QM atom coordinates
       _REAL_              :: qmvels(3,nqmatoms)   ! QM atom velocities (of previous step)
       integer, intent(in) :: nclatoms             ! Number of MM atoms
       _REAL_,  intent(in) :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
       _REAL_, intent(out) :: escf                 ! SCF energy
       _REAL_, intent(out) :: dxyzqm(3,nqmatoms)   ! SCF QM force
       _REAL_, intent(out) :: dxyzcl(3,nclatoms)   ! SCF MM force
       _REAL_              :: dipxyz(3), dipole    ! Dipole moment

       type(lio_nml_type), save     :: lio_nml
       logical, save                :: first_call = .true.
       integer                      :: nn,i
       integer                      :: printed =-1 ! Used to tell if we have printed this step yet 
                                                ! since the same step may be called multiple times
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! COMMANDS
!
! Replace with the following blocks between lines.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
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

!------------------------------------------------------------------------------!
! FFR - EHRENFEST ADAPTATION
!------------------------------------------------------------------------------!
    do nn=1,nqmatoms
       qmvels(1,nn)=x(lvel+3*qmmm_struct%iqmatoms(nn)-3)
       qmvels(2,nn)=x(lvel+3*qmmm_struct%iqmatoms(nn)-2)
       qmvels(3,nn)=x(lvel+3*qmmm_struct%iqmatoms(nn)-1)
    enddo

!    call SCF_in(escf,qmcoords,clcoords,nclatoms,dipxyz)
    call ehren_in( qmcoords, qmvels, clcoords, nclatoms, dipxyz, escf)
!------------------------------------------------------------------------------!
    escf=escf*CODATA08_AU_TO_KCAL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
