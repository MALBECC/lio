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
! Modified by Mariano C. González Lebrero - March 2013
! Modified by Francisco Ramirez           - May   2017
! Modified by Federico Pedron             - June  2020
!
! ----------------------------------------------------------------

  implicit none

  private
  public :: get_lio_forces

contains

  ! ---------------------------------
  ! Get QM energy and forces from LIO 
  ! ---------------------------------
  subroutine get_lio_forces(do_grad, nstep, ntpr_default, id, nqmatoms, &
                            qmcoords, nclatoms, clcoords, escf, dxyzqm, &
                            dxyzcl)

    use constants    , only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS, ZERO
    use file_io_dat  , only: mdin
    use memory_module, only: x, lvel
    use qmmm_module  , only: qmmm_struct, qmmm_nml

    implicit none
#   include "../include/md.h"

    integer, intent(in)  :: nqmatoms             ! Number of QM atoms
    _REAL_,  intent(in)  :: qmcoords(3,nqmatoms) ! QM atom coordinates
    integer, intent(in)  :: nclatoms             ! Number of MM atoms
    _REAL_,  intent(in)  :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
    _REAL_,  intent(out) :: escf                 ! SCF energy
    _REAL_,  intent(out) :: dxyzqm(3,nqmatoms)   ! SCF QM force
    _REAL_,  intent(out) :: dxyzcl(3,nclatoms)   ! SCF MM force
    _REAL_               :: dipxyz(3)            ! Dipole moment
    _REAL_               :: qmvels(3,nqmatoms)   ! QM atom velocities (of previous step)

    ! Unused dummy arguments.
    logical, intent(in) :: do_grad
    integer, intent(in) :: nstep
    integer, intent(in) :: ntpr_default
    character(len=3), intent(in) :: id

    logical, save :: first_call = .true.
    logical, save :: do_ehren_fsh = .false.
    integer       :: nn = 0

    ! Setup on first call.
    if ( first_call ) then
      first_call = .false.
      write (6,'(/,a,/)') 'Running QM/MM calculations with LIO'
      call init_lio_amber_new(nqmatoms, qmmm_struct%iqm_atomic_numbers, &
                              nclatoms, qmmm_nml%qmcharge, dt, mdin,    &
                              do_ehren_fsh)
    endif

    ! Get SCF density and energy.
    if (.not. do_ehren_fsh) then
      call SCF_in(escf, qmcoords, clcoords, nclatoms, dipxyz)      
    else
      do nn = 1, nqmatoms
        qmvels(1,nn) = x(lvel+3*qmmm_struct%iqmatoms(nn)-3)
        qmvels(2,nn) = x(lvel+3*qmmm_struct%iqmatoms(nn)-2)
        qmvels(3,nn) = x(lvel+3*qmmm_struct%iqmatoms(nn)-1)
      enddo
      call ehren_in( qmcoords, qmvels, clcoords, nclatoms, dipxyz, escf)
    endif
    escf = escf * CODATA08_AU_TO_KCAL

    ! Gradients/forces
    dxyzqm = 0.0D0
    dxyzcl = 0.0D0

    call dft_get_qm_forces(dxyzqm)
    call dft_get_mm_forces(dxyzcl,dxyzqm)

    dxyzqm = dxyzqm * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS
    dxyzcl = dxyzcl * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS

  end subroutine get_lio_forces
end module qm2_extern_lio_module
