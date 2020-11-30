#include "datatypes/datatypes.fh"
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% LIOMAIN.F90  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains the liomain subroutine, which performs several common     !
! procedures before and after calling SCF either by LIO alone or when          !
! performing in tantem with AMBER/GROMACS. Currently the routine:              !
! * Allocates matrices needed for SCF calculations.                            !
! * Calls SCF or SCFOP for closed/open shell calculations.                     !
! Other interfacial routines included are:                                     !
! * do_forces        (calculates forces/gradients)                             !
! * do_population_analysis (performs the analysis required)                    !
! * do_fukui()       (performs Fukui function calculation and printing)        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine liomain(E, dipxyz)
   use basis_data      , only: M, Nuc
   use cdft_data       , only: doing_cdft
   use cdft_subs       , only: cdft
   use cubegen         , only: cubegen_write, integrate_rho
   use cubegen_data    , only: cubegen_only
   use ehrensubs       , only: ehrendyn_main
   use fileio          , only: write_orbitals, write_orbitals_op
   use garcha_mod      , only: Smat, RealRho, OPEN, writeforces, energy_freq, &
                               NCO, restart_freq, npas, sqsm, natom,          &
                               hybrid_forces_props, doing_ehrenfest, Eorbs,   &
                               first_step, Eorbs_b, print_coeffs, MO_coef_at, &
                               MO_coef_at_b, Pmat_vec, nunp, r, d, Iz, pc,    &
                               rhoalpha, rhobeta
   use geometry_optim  , only: do_steep, doing_steep
   use mask_ecp        , only: ECP_init
   use tbdft_data      , only: MTB, tbdft_calc
   use tbdft_subs      , only: tbdft_init, tbdft_scf_output, write_rhofirstTB, &
                               tbdft_calibration
   use td_data         , only: tdrestart, timedep
   use time_dependent  , only: TD
   use typedef_operator, only: operator
   use dos_subs        , only: init_PDOS, build_PDOS, write_DOS, PDOS_finalize
   use excited_data    , only: excited_forces, pack_dens_exc
   use rhoint          , only: write1Drho
   use properties      , only: do_dipole, dipole
   use extern_functional_data, only: libint_inited

   implicit none
   LIODBLE  , intent(inout) :: E, dipxyz(3)

   type(operator) :: rho_aop, fock_aop, rho_bop, fock_bop
   integer        :: M_f, NCO_f, MM
   logical        :: calc_prop
   LIODBLE, allocatable :: Dens(:)

   call g2g_timer_sum_start("Total")
   npas = npas + 1

   ! TBDFT: Updating M and NCO for TBDFT calculations
   M_f   = M
   NCO_f = NCO
   if (tbdft_calc /= 0) then
      call tbdft_init(M, Nuc, OPEN)
      M_f   = M_f    + MTB
      NCO_f = NCO_f  + MTB / 2
   endif

   ! Libint initialization variable
   if ( libint_inited ) libint_inited = .false.

   if (.not.allocated(Smat))      allocate(Smat(M,M))
   if (.not.allocated(RealRho))   allocate(RealRho(M,M))
   if (.not.allocated(sqsm))      allocate(sqsm(M,M))
   if (.not.allocated(Eorbs))     allocate(Eorbs(M_f))
   if (.not.allocated(Eorbs_b))   allocate(Eorbs_b(M_f))

   call ECP_init()
   if (doing_steep()) then
      call do_steep(E)
   else if (doing_ehrenfest) then
      if (first_step) call SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
      call ehrendyn_main(E, dipxyz)
   else if (cubegen_only) then
      call cubegen_write(MO_coef_at)
   else if (tbdft_calc == 4) then
      call tbdft_calibration(E, fock_aop, rho_aop, fock_bop, rho_bop)
   else
      if (.not. tdrestart) then
         if (doing_cdft) then
            call CDFT(fock_aop, rho_aop, fock_bop, rho_bop, Pmat_vec, &
                      MO_coef_at, MO_coef_at_b, Smat, natom, M, NCO,  &
                      NCO+NUNP, OPEN)
         else
            call SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
         endif
      endif

      if (timedep == 1) then
         call TD(fock_aop, rho_aop, fock_bop, rho_bop)
      endif
   endif

   ! TBDFT calculations post SCF calculation:
   call tbdft_scf_output(OPEN)
   call write_rhofirstTB(M_f, OPEN)

   ! DOS and PDOS calculation post SCF calculation
   call init_PDOS(M_f)
   if (.not. OPEN) then
      call build_PDOS(MO_coef_at, Smat, M, M_f, Nuc)
      call write_DOS(M_f, Eorbs, 0)
   else
      call build_PDOS(MO_coef_at, Smat, M, M_f, Nuc)
      call write_DOS(M_f, Eorbs, 1)
      call build_PDOS(MO_coef_at_b, Smat, M, M_f, Nuc)
      call write_DOS(M_f, Eorbs_b,2)
   endif
   call PDOS_finalize()

   if ((restart_freq > 0) .and. (MOD(npas, restart_freq) == 0)) &
      call do_restart(88, Pmat_vec)

   ! Perform Mulliken and Lowdin analysis, get fukui functions and dipole.
   calc_prop = .false.
   if (((MOD(npas, energy_freq) == 0) .or. (hybrid_forces_props)) .and. &
       (.not. cubegen_only)) calc_prop = .true.

   ! Excited Properties
   MM = M*(M+1)/2
   allocate(Dens(MM))
   if ( excited_forces ) then
      Dens = pack_dens_exc
   else
      Dens = Pmat_vec
   endif

   if (calc_prop) then
      call cubegen_write(MO_coef_at(MTB+1:MTB+M,1:M))
      call do_population_analysis(Dens, rhoalpha, rhobeta)
      call do_fukui_calc()
      if (do_dipole()) call dipole(dipxyz, Pmat_vec, 2*NCO + nunp, &
                                   r, d, Iz, pc, .true., 1)
      if (writeforces) call do_forces(123)
      if (write1Drho) call integrate_rho()

      if (print_coeffs) then
         if (open) then
            call write_orbitals_op(M_f, NCO_f, NUnp, Eorbs, Eorbs_b, &
                                 MO_coef_at, MO_coef_at_b, 29)
         else
            call write_orbitals(M_f, NCO_f, Eorbs, MO_coef_at, 29)
         endif
      endif
   endif
   deallocate(Dens)

   call g2g_timer_sum_pause("Total")
end subroutine liomain

!%% DO_FORCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Calculates forces for QM and MM regions and writes them to output.           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_forces(uid)
    use garcha_mod, only: natom, nsol
    use fileio    , only: write_forces

    implicit none
    integer     , intent(in)  :: uid
    LIODBLE, allocatable :: dxyzqm(:,:), dxyzcl(:,:)

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

!%% DO_POPULATION_ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs the different population analyisis available.                       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_population_analysis(rho_tot, rho_a, rho_b)
   use garcha_mod      , only: Smat, Iz, sqsm, OPEN
   use properties      , only: print_mulliken, print_becke, print_lowdin, &
                               do_becke, do_mulliken, do_lowdin
   use basis_data      , only: M, Nuc, MM
   use ECP_mod         , only: ecpmode, IzECP
   use SCF_aux         , only: fix_densmat

   implicit none
   LIODBLE, intent(in)  :: rho_a(MM), rho_b(MM), rho_tot(MM)

   LIODBLE, allocatable :: rho_m(:,:), rho_mb(:,:)
   integer, allocatable :: true_iz(:)
   
   if ((.not. do_becke()) .and. (.not. do_mulliken()) &
                          .and. (.not. do_lowdin())) return

   call g2g_timer_sum_start('Population Analysis')
   
   ! Iz used to write the population file(s).
   allocate(true_iz(size(Iz,1)))
   true_iz = Iz
   if (ecpmode) true_iz = IzECP

   if ((do_lowdin()) .or. (do_mulliken())) then

      allocate(rho_m(M,M))
      allocate(rho_mb(1,1))

      if (open) then
         deallocate(rho_mb)
         allocate(rho_mb(M,M))

         call spunpack_rho('L', M, rho_a, rho_m)
         call spunpack_rho('L', M, rho_b, rho_mb)
      else 
         call spunpack_rho('L', M, rho_tot, rho_m)
      endif

      ! Performs Mulliken Population Analysis if required.
      if (do_mulliken()) then
         if (open) then
            call print_mulliken(rho_m, rho_mb, Smat, nuc, Iz, true_iz)
         else
            call print_mulliken(rho_m, Smat, nuc, Iz, true_iz)
         endif
      endif

      ! Performs Lowdin Population Analysis if required.
      if (do_lowdin()) then
         if (open) then
            call print_lowdin(rho_m, rho_mb, sqsm, nuc, Iz, true_iz)
         else
            call print_lowdin(rho_m, sqsm, nuc, Iz, true_iz)
         endif
      endif

      deallocate(rho_m, rho_mb)
   endif

   if (do_becke()) call print_becke(true_iz, open)

   deallocate(true_iz)
   call g2g_timer_sum_pause('Population Analysis')
endsubroutine do_population_analysis
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% do_fukui %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs Fukui function calls and printing.                                  !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_fukui_calc()
   use garcha_mod, only: MO_coef_at, MO_coef_at_b, NCO, Nunp, Smat, Eorbs, &
                         Eorbs_b, Iz, OPEN
   use properties, only: do_fukui, print_fukui
   use basis_data, only: Nuc
   use ECP_mod   , only: ecpmode, IzECP
   use tbdft_data, only: MTB, tbdft_calc
   
   implicit none
   integer :: NCO_f
   integer, allocatable :: true_iz(:)

   if (.not. do_fukui()) return
   NCO_f = NCO
   if (tbdft_calc /= 0) NCO_f = NCO_f + MTB/2

   allocate(true_iz(size(Iz,1)))
   true_iz = Iz
   if (ecpmode) true_iz = IzECP

   if (OPEN) then
      call print_fukui(MO_coef_at, MO_coef_at_b, NCO_f, NCO_f+NUNP, Nuc, &
                       Smat, Eorbs, Eorbs_b, Iz, true_iz)
   else
      call print_fukui(MO_coef_at, NCO_f, Nuc, Smat, Eorbs, Iz, true_iz)
   endif

   deallocate(true_iz)
    
end subroutine do_fukui_calc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% do_fukui() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs Fukui function calls and printing.                                  !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_restart(UID, rho_total)
   use garcha_mod , only: OPEN, NCO, NUNP, MO_coef_at, MO_coef_at_b, &
                          rhoalpha, rhobeta
   use basis_data , only: M, MM, indexii
   use fileio_data, only: rst_dens
   use fileio     , only: write_coef_restart, write_rho_restart
   use tbdft_data,  only: MTB, tbdft_calc

   implicit none
   integer         , intent(in) :: UID
   LIODBLE, intent(in) :: rho_total(MM)
   LIODBLE, allocatable :: coef(:,:), coef_b(:,:), tmp_rho(:,:), &
                                    tmp_rho_b(:,:)
   integer :: NCOb, icount, jcount
   integer :: NCO_f, i0
!TBDFT: Updating M for TBDFT calculations
   if (tbdft_calc/=0) then
      NCO_f = NCO + MTB/2
      i0    = MTB
   else
      NCO_f = NCO
      i0    = 0
   end if

   if ( rst_dens == 2 ) then
      allocate(tmp_rho(M,M))
      if (.not. OPEN) then
         call spunpack('L', M, rho_total, tmp_rho)
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
      allocate(coef(M, NCO_f))
      do icount=1, M
      do jcount=1, NCO_f
         coef(indexii(icount), jcount) = MO_coef_at(i0+icount,jcount)
      enddo
      enddo
      if (OPEN) then
         NCOb = NCO_f + NUNP
         allocate(coef_b(M, NCOb))

         do icount=1, M
         do jcount=1, NCOb
            coef_b(indexii(icount), jcount) = MO_coef_at_b(i0+icount,jcount)
         enddo
         enddo
         call write_coef_restart(coef, coef_b, M, NCO_f, NCOb, UID)
         deallocate(coef_b)
      else
         call write_coef_restart(coef, M, NCO_f, UID)
      endif
      deallocate(coef)
   endif

   return
end subroutine do_restart
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
