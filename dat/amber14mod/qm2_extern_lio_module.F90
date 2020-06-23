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
    use memory_module, only: x, lvel, i04
    use parms        , only: cn1, cn2, nttyp
    use qmmm_module  , only: qmmm_struct, qmmm_nml

    use memory_module, only: glob_int_array => ix
    ! We only need this for LJS, ix should be explored from I04 from memory.h.

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
    logical, save :: do_ljswitch  = .false. ! True if doing LJ Switch corrections.
    integer       :: nn = 0

    integer :: tmp_idx ! Temporary index

    ! Variables for LJ Switch (if used)
    integer :: ntyp                  ! Number of LJ types.
    _REAL_ , allocatable :: sig(:)   ! List of sigmas given to LIO
    _REAL_ , allocatable :: eps(:)   ! List of epsilons given to LIO
    integer, allocatable :: qm_types(:)
    integer, allocatable :: mm_types(:)

    ! Setup on first call.
    if ( first_call ) then
      first_call = .false.
      write (6,'(/,a,/)') 'Running QM/MM calculations with LIO'
      call init_lio_amber_new(nqmatoms, qmmm_struct%iqm_atomic_numbers, &
                              nclatoms, qmmm_nml%qmcharge, dt, mdin,    &
                              do_ehren_fsh, do_ljswitch)

       if (do_ljswitch) then
        ! Here we recover the number of types, since 
        ! nttyp is ntyp * (ntyp + 1) / 2
        ntyp = floor(sqrt(2.0D0 * nttyp))

        ! Gets sigma and epsilon from CN1 (A) and CN2 (B).
        ! Epsilon is stored as 4e.
        allocate(sig(ntyp), eps(ntyp))
        sig = 0.0D0; eps = 0.0D0
        do nn = 1, ntyp
          tmp_idx = nn * (nn + 1) / 2
          if ((abs(cn1(tmp_idx)) > 1D-5) .and. (abs(cn2(tmp_idx)) > 1D-5)) then
            eps(nn) = cn2(tmp_idx) * cn2(tmp_idx) / cn1(tmp_idx)
            sig(nn) = ( cn1(tmp_idx) / cn2(tmp_idx) ) ** (1.0D0 / 6.0D0)
          endif
        enddo

        eps = eps / CODATA08_AU_TO_KCAL
        sig = sig * CODATA08_A_TO_BOHRS
        call ljs_set_params(ntyp, eps, sig)

        deallocate(sig, eps)
      endif
    endif

    ! Get SCF density and energy.
    if (.not. do_ehren_fsh) then

      if (do_ljswitch) then
        allocate(qm_types(nqmatoms), mm_types(nclatoms))
        do nn = 1, nqmatoms
          tmp_idx = I04 + qmmm_struct%iqmatoms(nn) -1
          qm_types(nn) = glob_int_array(tmp_idx)
        enddo

        do nn = 1, nclatoms
          tmp_idx = I04 + qmmm_struct%qm_mm_pair_list(nn) -1
          mm_types(nn) = glob_int_array(tmp_idx)
        enddo
        
        call ljs_settle_mm(qm_types, mm_types, qmcoords * CODATA08_A_TO_BOHRS,&
                           clcoords * CODATA08_A_TO_BOHRS, nqmatoms, nclatoms)
        deallocate(qm_types, mm_types)
      endif

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

    if (do_ljswitch) then
      call ljs_substract_mm(escf, dxyzqm, dxyzcl, qmcoords*CODATA08_A_TO_BOHRS,&
                            clcoords * CODATA08_A_TO_BOHRS, nqmatoms, nclatoms)
    endif

    dxyzqm = dxyzqm * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS
    dxyzcl = dxyzcl * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS

  end subroutine get_lio_forces
end module qm2_extern_lio_module
