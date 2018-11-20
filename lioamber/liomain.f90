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
    use garcha_mod, only: M, Smat, RealRho, OPEN, writeforces, energy_freq,   &
                          restart_freq, npas, sqsm, mulliken, lowdin, dipole, &
                          doing_ehrenfest, first_step, Eorbs, Eorbs_b, fukui, &
                          print_coeffs, steep,       MO_coef_at, MO_coef_at_b,&
                          NUnp, NCO
    use ecp_mod   , only: ecpmode, IzECP
    use ehrensubs , only: ehrendyn_main
    use fileio    , only: write_orbitals, write_orbitals_op

    implicit none
    REAL*8, intent(inout) :: dipxyz(3), E

    call g2g_timer_sum_start("Total")

    if (.not.allocated(Smat))      allocate(Smat(M,M))
    if (.not.allocated(RealRho))   allocate(RealRho(M,M))
    if (.not.allocated(sqsm))      allocate(sqsm(M,M))
    if (.not.allocated(Eorbs))     allocate(Eorbs(M))
    if (.not.allocated(Eorbs_b))   allocate(Eorbs_b(M))

    if (steep) then
      call do_steep(E)
    end if

    if ( doing_ehrenfest ) then
       if ( first_step ) call SCF( E, dipxyz )
       call ehrendyn_main( E, dipxyz )
    else
       call SCF(E)
    endif

    if ( (restart_freq.gt.0) .and. (MOD(npas, restart_freq).eq.0) ) &
       call do_restart(88)

    ! Perform Mulliken and Lowdin analysis, get fukui functions and dipole.
    if (MOD(npas, energy_freq).eq.0) then
        if (mulliken.or.lowdin) call do_population_analysis()
        if (dipole) call do_dipole(dipxyz, 69)
        if (fukui) call do_fukui()

        if (writeforces) then
            if (ecpmode) stop "ECP does not feature forces calculation."
            call do_forces(123)
        endif
        if (print_coeffs) then
           if (open) then
             call write_orbitals_op(M, NCO, NUnp, Eorbs, Eorbs_b, MO_coef_at,  &
                                    MO_coef_at_b, 29)
          else
             call write_orbitals(M, NCO, Eorbs, MO_coef_at, 29)
          endif
        endif
    endif

    call g2g_timer_sum_pause("Total")

    return
end subroutine liomain

!%% DO_FORCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Calculates forces for QM and MM regions and writes them to output.           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_forces(uid)
    use garcha_mod, only: natom, nsol
    use fileio    , only: write_forces

    implicit none
    integer, intent(in) :: uid
    integer             :: k
    real*8, allocatable :: dxyzqm(:,:), dxyzcl(:,:)

    call g2g_timer_start('Forces')
    open(unit=uid, file='forces')

    allocate ( dxyzqm(3, natom) )
    dxyzqm = 0.0D0

    call dft_get_qm_forces(dxyzqm)
    if (nsol.gt.0) then
        allocate ( dxyzcl(3, natom+nsol) )
        dxyzcl = 0.0D0
        call dft_get_mm_forces(dxyzcl, dxyzqm)
    endif

    call write_forces(dxyzqm, natom, 0, uid)
    deallocate (dxyzqm)

    if(nsol.gt.0) then
        call write_forces(dxyzcl, nsol, 0, uid)
        deallocate (dxyzcl)
    endif
    call g2g_timer_stop('Forces')

    return
end subroutine do_forces
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% DO_DIPOLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Sets variables up and calls dipole calculation.                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_dipole(dipxyz, uid)
    use fileio, only: write_dipole
    implicit none
    integer, intent(in)    :: uid
    real*8 , intent(inout) :: dipxyz(3)
    real*8                 :: u

    call g2g_timer_start('Dipole')
    call dip(dipxyz)
    u = sqrt(dipxyz(1)**2 + dipxyz(2)**2 + dipxyz(3)**2)

    call write_dipole(dipxyz, u, uid, .true.)
    call write_dipole(dipxyz, u, uid, .false.)
    call g2g_timer_stop('Dipole')

    return
end subroutine do_dipole
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% DO_POPULATION_ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs the different population analyisis available.                       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_population_analysis()
   use garcha_mod, only: RMM, Smat, RealRho, M, Enucl, Nuc, Iz, natom, &
                         mulliken, lowdin, sqsm, d, r, Md, ntatom
   use ECP_mod   , only: ecpmode, IzECP
   use faint_cpu , only: int1
   use SCF_aux   , only: fix_densmat
   use fileio    , only: write_population

   implicit none
   integer :: MM, M11, M5, IzUsed(natom), kk
   real*8  :: q(natom), En

   ! Needed until we dispose of RMM.
   MM = M*(M+1)/2 ; M5 = 1 + MM*2; M11 = M5+MM+Md*(Md+1)

   ! Decompresses and fixes S and RealRho matrixes, which are needed for
   ! population analysis.
   call int1(En, RMM(M5:M5+MM), RMM(M11:M11+MM), Smat, d, r, Iz, natom, ntatom)
   call spunpack('L',M,RMM(M5),Smat)
   call spunpack('L',M,RMM(1),RealRho)
   call fix_densmat(RealRho)

   ! Initial nuclear charge for Mulliken
   do kk=1,natom
       q(kk) = real(Iz(kk))
   enddo

   ! Iz used to write the population file.
   IzUsed = Iz
   if (ecpmode) IzUsed = IzECP

   ! Performs Mulliken Population Analysis if required.
   if (mulliken) then
       call g2g_timer_start('Mulliken')
       call mulliken_calc(natom, M, RealRho, Smat, Nuc, q)
       call write_population(natom, IzUsed, q, 0, 85)
       call g2g_timer_stop('Mulliken')
   endif

   ! Performs Löwdin Population Analysis if required.
   if (lowdin) then
       call g2g_timer_start('Lowdin')
       call lowdin_calc(M, natom, RealRho, sqsm, Nuc, q)
       call write_population(natom, IzUsed, q, 1, 85)
       call g2g_timer_stop('Lowdin')
   endif

   return
endsubroutine do_population_analysis
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% DO_FUKUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs Fukui function calls and printing.                                  !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_fukui()
    use garcha_mod, only: MO_coef_at, MO_coef_at_b, NCO, Nunp, M, natom, Nuc, &
                          Smat, Eorbs, Eorbs_b, Iz, OPEN
    use fileio    , only: write_fukui
    implicit none
    double precision              :: softness
    double precision, allocatable :: fukuim(:), fukuin(:), fukuip(:)

    allocate(fukuim(natom), fukuin(natom), fukuip(natom))

    call g2g_timer_start("Fukui")
    if (OPEN) then
        call fukui_calc_os(MO_coef_at, MO_coef_at_b, NCO, NCO+Nunp, M, natom, &
                           Nuc, Smat, fukuim, fukuip, fukuin, Eorbs, Eorbs_b)
        call get_softness(Eorbs(NCO-1), Eorbs(NCO), Eorbs_b(NCO+NUNP-1), &
                          Eorbs(NCO+NUNP), softness)
        call write_fukui(fukuim, fukuip, fukuin, natom, Iz, softness)
    else
        call fukui_calc(MO_coef_at, NCO, M, natom, Nuc, Smat, fukuim, fukuip, &
                        fukuin, Eorbs)
        call get_softness(Eorbs(NCO-1), Eorbs(NCO), Eorbs(NCO-1), Eorbs(NCO), &
                          softness)
        call write_fukui(fukuim, fukuip, fukuin, natom, Iz, softness)
    endif
    call g2g_timer_stop("Fukui")

    deallocate(fukuim, fukuin, fukuip)
end subroutine do_fukui
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% DO_FUKUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs Fukui function calls and printing.                                  !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_restart(UID)
   use garcha_mod , only: RMM, OPEN, NCO, NUNP, M, MO_coef_at, MO_coef_at_b, &
                          indexii, rhoalpha, rhobeta
   use fileio_data, only: rst_dens
   use fileio     , only: write_coef_restart, write_rho_restart
   implicit none
   integer, intent(in) :: UID
   integer             :: NCOb, icount, jcount, coef_ind
   real*8, allocatable :: coef(:,:), coef_b(:,:), tmp_rho(:,:), tmp_rho_b(:,:)

   if ( rst_dens .eq. 2 ) then
      allocate(tmp_rho(M,M))
      if (.not. OPEN) then
         call spunpack('L', M, RMM(1), tmp_rho)
         call write_rho_restart(tmp_rho, M, uid)
      else
         allocate(tmp_rho_b(M,M))
         call spunpack('L', M, rhoalpha, tmp_rho)
         call spunpack('L', M, rhobeta , tmp_rho_b)
         call write_rho_restart(tmp_rho, tmp_rho_b, M, uid)
         deallocate(tmp_rho_b)
      endif
      deallocate(tmp_rho)
   else
      allocate(coef(M, NCO))
      do icount=1, M
      do jcount=1, NCO
         coef_ind = icount + M*(jcount-1)
         coef(indexii(icount), jcount) = MO_coef_at(coef_ind)
      enddo
      enddo

      if (OPEN) then
         NCOb = NCO + NUNP
         allocate(coef_b(M, NCOb))

         do icount=1, M
         do jcount=1, NCOb
            coef_ind = icount + M*(jcount-1)
            coef_b(indexii(icount), jcount) = MO_coef_at_b(coef_ind)
         enddo
         enddo

         call write_coef_restart(coef, coef_b, M, NCO, NCOb, UID)
         deallocate(coef_b)
      else
         call write_coef_restart(coef, M, NCO, UID)
      endif
      deallocate(coef)
   endif

   return
end subroutine do_restart
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
