!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% LIOMAIN.F90  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains the liomain subroutine, which performs several common     !
! procedures before and after calling SCF either by LIO alone or when          !
! performing in tantem with AMBER/GROMACS. Currently the routine:              !
! * Allocates matrices needed for SCF calculations.                            !
! * Calls SCF or SCFOP for closed/open shell calculations.                     !
! Other interfacial routines included are:                                     !
! * do_forces        (calculates forces/gradients)                             !
! * do_dipole        (calculates dipole moment)                                !
! * do_population_analysis (performs the analysis required)                    !
! * do_fukui         (performs Fukui function calculation and printing)        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
 
subroutine liomain(E, dipxyz)
    use garcha_mod, only : M, Smat, RealRho, OPEN, writeforces, energy_freq,   &
                           restart_freq, npas, sqsm, mulliken, lowdin, dipole, &
                           do_ehrenfest, first_step,                           &
                           Eorbs, fukui, print_coeffs
    use ecp_mod   , only : ecpmode, IzECP
    use ehrenfest,  only : ehrendyn
 
    implicit none
    REAL*8, intent(inout) :: dipxyz(3), E

    if (.not.allocated(Smat))    allocate(Smat(M,M))
    if (.not.allocated(RealRho)) allocate(RealRho(M,M))
    if (.not.allocated(sqsm))    allocate(sqsm(M,M))
    if (.not.allocated(Eorbs))   allocate(Eorbs(M))

!------------------------------------------------------------------------------!
! FFR - Option to do ehrenfest
    if ( do_ehrenfest ) then
       if ( first_step ) call SCF( E, dipxyz )
       call ehrendyn( E, dipxyz )

    else
       if(OPEN) then
          if (ecpmode) stop "ECP is unavailable for Open Shell systems."
          call SCFOP(E, dipxyz)
       else
          call SCF(E)
       endif

    endif

    ! Perform Mulliken and Lowdin analysis, get fukui functions and dipole.
    if (MOD(npas, energy_freq).eq.0) then

        if (mulliken.or.lowdin) call do_population_analysis()

        if (dipole) call do_dipole(dipxyz, 69)
  
        if (fukui) call do_fukui()

        if (writeforces) then
            if (ecpmode) stop "ECP does not feature forces calculation."
            call do_forces(123)
        endif

        if (print_coeffs) call write_orbitals(29)
    endif

    if ((restart_freq.gt.0).and.(MOD(npas, restart_freq).eq.0)) call write_restart(88)

    return
end subroutine liomain


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DO_FORCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Calculates forces for QM and MM regions and writes them to output.           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_forces(uid)

    use garcha_mod, only : natom, nsol

    implicit none
    integer, intent(in) :: uid
    integer             :: k
    real*8, allocatable :: dxyzqm(:,:), dxyzcl(:,:)

    call g2g_timer_start('Forces')
    open(unit=uid, file='forces')

    allocate ( dxyzqm(3, natom) )
    dxyzqm = 0.0

    call dft_get_qm_forces(dxyzqm)
    if (nsol.gt.0) then
        allocate ( dxyzcl(3, natom+nsol) )
        dxyzcl = 0.0
        call dft_get_mm_forces(dxyzcl, dxyzqm)
    endif

    call write_forces(dxyzqm, natom, 0, uid)
    deallocate (dxyzqm)

    if(nsol.gt.0) then
        call write_forces(dxyzcl, nsol, natom, uid)
        deallocate (dxyzcl)
    endif
    call g2g_timer_stop('Forces')

    return
end subroutine do_forces
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DO_DIPOLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Sets variables up and calls dipole calculation.                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_dipole(dipxyz, uid)
    implicit none
    integer, intent(in)    :: uid
    real*8 , intent(inout) :: dipxyz(3)
    real*8                 :: u

    call g2g_timer_start('Dipole')
    call dip(dipxyz)
    u = sqrt(dipxyz(1)**2 + dipxyz(2)**2 + dipxyz(3)**2)

    call write_dipole(dipxyz, u, uid)
    call g2g_timer_stop('Dipole')

    return
end subroutine do_dipole
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DO_POPULATION_ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs the different population analyisis available.                       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_population_analysis()
   use garcha_mod, only : RMM, Smat, RealRho, M, Enucl, Nuc, Iz, natom, &
                          mulliken, lowdin, sqsm
   use ECP_mod   , only : ecpmode, IzECP

   implicit none
   integer :: M1, M5, IzUsed(natom), kk
   real*8  :: q(natom)

   ! Needed until we dispose of RMM.
   M1=1 ; M5=1+M*(M+1)

   ! Iz used to write the population file.
   IzUsed = Iz
   if (ecpmode) IzUsed = IzECP

   ! Decompresses and fixes S and RealRho matrixes, which are needed for
   ! population analysis.
   call int1(Enucl)
   call spunpack('L',M,RMM(M5),Smat)
   call spunpack('L',M,RMM(M1),RealRho)
   call fixrho(M,RealRho)

   do kk=1,natom
       q(kk) = real(IzUsed(kk))
   enddo

   ! Performs Mulliken Population Analysis if required.
   if (mulliken) then
       call g2g_timer_start('Mulliken')
       call mulliken_calc(natom, M, RealRho, Smat, Nuc, q)
       call write_population(85, natom, Iz, q, 0)
       call g2g_timer_stop('Mulliken')
   endif

   ! Performs Löwdin Population Analysis if required.
   if (lowdin) then
       call g2g_timer_start('Lowdin')
       call lowdin_calc(M, natom, RealRho, sqsm, Nuc, q)
       call write_population(85, natom, Iz, q, 1)
       call g2g_timer_stop('Lowdin')
   endif

   return
endsubroutine do_population_analysis
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DO_FUKUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs Fukui function calls and printing.                                  !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_fukui()
    use garcha_mod, only : X, NCO, M, natom, Nuc, Smat, Eorbs, Iz, OPEN

    implicit none
    real*8  :: fukuim(natom), fukuin(natom), fukuip(natom), softness


    call g2g_timer_start("Fukui")
    if (OPEN) then
        ! TO-DO
        ! Add call fukui_calc_OS when openshell is available.
    else
        call fukui_calc(X(1,M*2+1), NCO, M, natom, Nuc, Smat, fukuim, fukuip, &
                        fukuin, Eorbs)
        call get_softness(Eorbs(NCO-1), Eorbs(NCO), Eorbs(NCO-1), Eorbs(NCO), &
                          softness)
        call write_fukui(fukuim, fukuip, fukuin, natom, Iz, softness)
    endif
    call g2g_timer_stop("Fukui")

    return
end subroutine do_fukui
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

