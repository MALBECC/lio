#include "datatypes/datatypes.fh"
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% REAL TIME-TDDFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! 1992 - Dario Estrin (source code)                                            !
! 2012 - Nano, Dario, Uriel, Damian (main routine, external fields)            !
! 2015 - Fran (magnus)                                                         !
! 2017 - Charly, Gonzalo, Fede (reformating)                                   !
!                                                                              !
! This subroutine takes the converged density matrix from an SCF calculation   !
! and evolves it in time. In the input file the total number of propagation    !
! steps is specified (nstep) as well as the time of each evolution step        !
! (tdstep).                                                                    !
! This implementation has two alternatives to evolve the density in time. The  !
! first one (propagator=1) is the Verlet algorithm that uses a convination of  !
! Liouville von Newmann expresion for the time derivative of the density matrix!
! and a first order Taylor expansion of the density matrix. The second one     !
! (propagator=2) is the Magnus propagation scheme that uses Backer Campbell    !
! Hausdorff (BCH) formula. For this reason when Magnus is used the number of   !
! total conmutators in the BCH espansion has to be specified (NBCH,            !
! default=10).                                                                 !
! A narrow gaussian type electric field can be introduced during the time      !
! evolution in order to excite all electronic frequencies with the same        !
! intensity. Once this perturbation is turned on (Field=t, exter=t) each       !
! component of the external electric field has to be specified in the input    !
! file (Fx,Fy,Fz).                                                             !
! In each step of the propagation the cartesian components of the sistems      !
! dipole are stored in files x.dip, y.dip, z.dip.                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

module td_data
   implicit none
   integer :: td_rst_freq = 500
   integer :: timedep    = 0
   integer :: ntdstep    = 0
   integer :: td_do_pop  = 0
   integer :: td_do_opop = 0
   integer :: td_eu_step = 0      !Include Euler steps each td_eu_step steps
                                  !only in Verlet propagator.
   LIODBLE :: tdstep     = 2.0D-3
   logical :: tdrestart  = .false.
   logical :: writedens  = .false.
   LIODBLE :: pert_time  = 2.0D-1

   integer :: td_rho_purify = 0
   logical :: td_check_rho  = .true. ! Verifies rho idempotency.
end module td_data

module time_dependent
   implicit none
contains

subroutine TD(fock_aop, rho_aop, fock_bop, rho_bop)
   use garcha_mod    , only: NBCH, propagator, Iz, igrid2, r, nsol, pc, Smat, &
                             MEMO, ntatom, sqsm, OPEN, natom, d, rhoalpha,    &
                             rhobeta, Fmat_vec, Fmat_vec2, Ginv_vec, Hmat_vec,&
                             Gmat_vec, Pmat_vec, NCO, nunp
   use basis_data    , only: M, Nuc, MM
   use basis_subs    , only: neighbour_list_2e
   use td_data       , only: td_rst_freq, tdstep, ntdstep, tdrestart, &
                             writedens, pert_time, td_eu_step, td_do_opop, &
                             td_do_pop, td_rho_purify, td_check_rho
   use field_data    , only: field, fx, fy, fz
   use field_subs    , only: field_setup_old, field_finalize
   use transport_data, only: transport_calc
   use transport_subs, only: transport_rho_trace, transport_generate_rho,      &
                             transport_init, transport_population
   use tbdft_data     , only: tbdft_calc, MTBDFT, MTB
   use tbdft_subs     , only: tbdft_td_init, tbdft_td_output
   use fileio         , only: write_td_restart_verlet, write_td_restart_magnus, &
                              read_td_restart_verlet , read_td_restart_magnus,  &
                              write_energies, movieprint
   use fileio_data     , only: verbose
   use typedef_operator, only: operator
   use typedef_cumat   , only: cumat_x, cumat_r
   use faint_cpu       , only: int2
   use ceed_subs       , only: ceed_init

   implicit none
!carlos: Operator inserted for TD
   type(operator), intent(inout)           :: rho_aop, fock_aop
   type(operator), intent(inout), optional :: rho_bop, fock_bop

   LIODBLE  :: E, En, E1, E2, E1s, Eexact, Es, Ens = 0.0D0, Ex, t, dt_magnus, dt_lpfrg
   integer :: M2, igpu, istep
   integer :: lpfrg_steps = 200, chkpntF1a = 185, chkpntF1b = 195
   logical :: is_lpfrg = .false. , fock_restart = .false.
   character(len=20) :: restart_filename

   LIODBLE , allocatable, dimension(:)   :: factorial
   LIODBLE , allocatable, dimension(:,:) :: overlap, Smat_initial
!carlos: the next variables have 3 dimensions, the 3th one is asociated with the
!        spin number. This one will have the value of 1 for Close Shell and 2
!        for Open shell. Spin=1 will always be reffered to alpha and spin=2,
!        beta.

   LIODBLE , allocatable, dimension(:,:,:) :: fock, F1a, F1b
   LIODBLE , allocatable, dimension(:,:,:) :: fock_0

! Precision options.
   TDCOMPLEX :: Im
   TDCOMPLEX :: liocmplx
   TDCOMPLEX, allocatable, dimension(:,:,:) :: rho, rho_aux, rhonew, rhold
   TDCOMPLEX, allocatable, dimension(:,:,:) :: rho_0

   type(cumat_r) :: Xmat
   type(cumat_x) :: Xtrans, Ymat

!TBDFT: M_f controls de size of the bigest matrices for TBDFT, ii and jj are only
!counters, and traza is for the control of the trace of density matrix
   integer :: M_f
!carlos: Open Shell variables

   integer :: dim3
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% TD INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   call g2g_timer_start('TD')
   call g2g_timer_start('td-inicio')
   open(unit = 134, file = "dipole_moment_td")

!------------------------------------------------------------------------------!
!TBDFT: defining TBDFT matrix size

      if (tbdft_calc /= 0) then
         M_f = MTBDFT
      else
         M_f = M
      end if
!------------------------------------------------------------------------------!
!carlos: dim3 is the 3th dimension of all matrix involve in open shell
!        calculation

   M2  = 2*M
   Im = liocmplx(0.0D0,2.0D0)
   dim3 = 1
   if (open) dim3 = 2

!------------------------------------------------------------------------------!

   ! Checks and performs allocations.
   call td_allocate_all(M_f, M, dim3, NBCH, propagator, F1a, F1b, fock, rho,   &
                        rho_aux, rhold, rhonew, rho_0, fock_0, sqsm, factorial,&
                        Smat_initial)

   ! Initialises propagator-related parameters and other variables.
   call td_initialise(propagator, tdstep, NBCH, dt_lpfrg, dt_magnus, factorial)
!TBDFT:TD restart is not working with TBDFT
   ! TD restart reading.
   if (tdrestart) then
      if (OPEN) then
         restart_filename = 'td_a_in.restart'
         if (propagator.eq.2) then
            call read_td_restart_magnus(rho_0(:,:,1), F1a(:,:,1), F1b(:,:,1),  &
                                        M_f, restart_filename, fock_restart)
         else
            call read_td_restart_verlet(rho_0(:,:,1), M_f, restart_filename)
         endif
         restart_filename = 'td_b_in.restart'
         if (propagator.eq.2) then
            call read_td_restart_magnus(rho_0(:,:,2), F1a(:,:,2), F1b(:,:,2),  &
                                        M_f, restart_filename, fock_restart)
         else
            call read_td_restart_verlet(rho_0(:,:,2), M_f, restart_filename)
         endif
         call sprepack_ctr('L', M, rhoalpha, rho_0(MTB+1:MTB+M,MTB+1:MTB+M,1))
         call sprepack_ctr('L', M, rhobeta , rho_0(MTB+1:MTB+M,MTB+1:MTB+M,2))
         Pmat_vec = rhoalpha + rhobeta
      else
         restart_filename = 'td_in.restart'
         if (propagator.eq.2) then
            call read_td_restart_magnus(rho_0(:,:,1), F1a(:,:,1), F1b(:,:,1),  &
                                        M_f, restart_filename, fock_restart)
         else
            call read_td_restart_verlet(rho_0(:,:,1), M_f, restart_filename)
         endif
         call sprepack_ctr('L', M, Pmat_vec, rho_0(MTB+1:MTB+M,MTB+1:MTB+M,1))
      end if
   else
      ! Read the density matrix stored in Pmat_vec(1,2,3,...,MM) into rho matrix.
      !For Open Shell, this is read from rhoalpha and rhobeta
      if (OPEN) then
         call spunpack_rtc('L', M, rhoalpha, rho_0(:,:,1))
         call spunpack_rtc('L', M, rhobeta, rho_0(:,:,2))
      else
         call spunpack_rtc('L', M, Pmat_vec, rho_0(:,:,1))
       end if
   endif

!------------------------------------------------------------------------------!
! TBDFT: Initialize TBDFT variables for TD
   if (tbdft_calc /= 0) then
      if(OPEN) then
         call tbdft_td_init (M , rho, rho_0, 2)
      else
         call tbdft_td_init (M , rho, rho_0, 1)
      end if
   else
      rho=rho_0
   end if
!carlos: storing rho AO data in Operator

   call rho_aop%Sets_dataC_AO(rho(:,:,1))
   if (OPEN) call rho_bop%Sets_dataC_AO(rho(:,:,2))
!------------------------------------------------------------------------------!

   ! Proper TD calculation start.
   if (verbose .gt. 2) write(*,'(A)') 'Starting TD calculation.'
   ! Create integration grid for XC, assigning points to groups (spheres/cubes)
   ! and significant functions to groups, also calculating point weights.
   if (field) call field_setup_old(pert_time, 1, fx, fy, fz)
   call neighbour_list_2e(natom, ntatom, r, d)
   call td_integration_setup(igrid2, igpu, Iz)
   call td_integral_1e(E1, En, E1s, Ens, MM, igpu, nsol, Pmat_vec, Fmat_vec, &
                       Hmat_vec, r, pc, ntatom, natom, Smat, d, Iz)
   call spunpack('L', M, Fmat_vec, Smat_initial)

   ! Initialises transport if required.
   if (transport_calc) call transport_init(M, dim3, natom, Nuc, Fmat_vec, &
                                           overlap, rho,OPEN)

   ! Diagonalizes Smat and calculates the base change matrices (x,y,Xtrans)
   call td_overlap_diag(M_f, M, Smat, Xmat, Xtrans, Ymat)

   call ceed_init(M, OPEN, r, d, natom, ntatom, propagator)

   call rho_aop%BChange_AOtoON(Ymat, M_f)
   if (OPEN) call rho_bop%BChange_AOtoON(Ymat, M_f)

   if (td_rho_purify > 0) then
      call rho_aop%purify_ON(OPEN)
      if (OPEN) call rho_bop%purify_ON(OPEN)
   endif
   if (td_check_rho) then
      call rho_aop%check_idempotency_ON(OPEN)
      if (OPEN) call rho_bop%check_idempotency_ON(OPEN)
   endif
   
   ! Precalculate three-index (two in MO basis, one in density basis) matrix
   ! used in density fitting /Coulomb F element calculation here (t_i in Dunlap)
   call int2(Gmat_vec, Ginv_vec, r, d, ntatom)
   call td_coulomb_precalc(igpu, MEMO, r, d, natom, ntatom)

   ! Recalculate maximum number of TD steps.
   if ((propagator .eq. 2) .and. (.not. fock_restart)) then
      if (verbose .gt. 2) write(*,'(A)') "  Magnus propagator selected. &
                                         &Adjusting maximum number of TD steps."
      ntdstep = ntdstep + 9 * lpfrg_steps / 10
   endif


   call g2g_timer_stop('td-inicio')
   ! End of TD initialization.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% TD EVOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   do 999 istep = 1, ntdstep
      call g2g_timer_start('TD step')
      call g2g_timer_sum_start("TD - TD Step")
      ! Checks if step is a leapfrog step.
      call td_check_prop(is_lpfrg, propagator, istep, lpfrg_steps, &
                         fock_restart, verbose)
      call td_get_time(t, tdstep, istep, propagator, is_lpfrg)

      call g2g_timer_sum_start("TD - TD Step Energy")

      call td_calc_energy(E, E1, E2, En, Ex, Es, Eexact, MM, Pmat_vec, Fmat_vec,&
                          Fmat_vec2, Gmat_vec, Ginv_vec, Hmat_vec, is_lpfrg,    &
                          t / 0.024190D0, M, open, r, d, natom, ntatom, MEMO,   &
                          2*NCO+NUNP, Iz, pc)

      call g2g_timer_sum_pause("TD - TD Step Energy")
      if (verbose > 2) write(*,'(A,I6,A,F12.6,A,F12.6,A)') "  TD Step: ", &
                                istep, " - Time: ", t, " fs - Energy : ", E, &
                                " A.U."

      ! Prints orbital population if required.
      call td_orbital_population(rho_aop, rho_bop, OPEN, istep, td_do_opop, &
                                 propagator, is_lpfrg, M)

      ! Verlet or Magnus Propagation
      ! In Verlet, we obtain the Fock matrix in the molecular orbital (MO)
      ! basis, where U matrix with eigenvectors of S, and s is vector with
      ! eigenvalues. In the first step of the propagation we extrapolate rho
      ! back in time using Verlet algorithm to calculate rhold.
      ! Both propagations include transport.

      ! After propagation, we transform the density to the atomic orbital basis
      ! and take the real part of it. The imaginary part of the density can be
      ! discarded since for a basis set of purely real functions the fock matrix
      ! is real and symetric, and depends only on the real part of the complex
      ! density matrix. (This will not be true in the case of hybrid
      ! functionals)

      call g2g_timer_sum_start("TD - Propagation")

      if (is_lpfrg) then
         call td_bc_fock(M_f, M, Fmat_vec, fock_aop, Xmat, istep, &
                         t/0.024190D0)
         if (OPEN) then
            call td_bc_fock(M_f, M, Fmat_vec2, fock_bop, Xmat, istep, &
                            t/0.024190D0)
            call td_verlet(M, M_f, dim3, OPEN, fock_aop, rhold, rho_aop, &
                           rhonew, istep, td_eu_step, Im, dt_lpfrg,            &
                           transport_calc, natom, Nuc, Iz, overlap, sqsm, Ymat,&
                           Xtrans, fock_bop, rho_bop)
         else
            call td_verlet(M, M_f, dim3, OPEN, fock_aop, rhold, rho_aop, &
                           rhonew, istep, td_eu_step, Im, dt_lpfrg,            &
                           transport_calc, natom, Nuc, Iz, overlap, sqsm, Ymat,&
                           Xtrans)
         end if

         if (propagator == 2) then
            call fock_aop%Gets_data_ON(fock(:,:,1))
            if (OPEN) call fock_bop%Gets_data_ON(fock(:,:,2))
            if (istep == chkpntF1a) F1a = fock
            if (istep == chkpntF1b) F1b = fock
         endif
      else
         if(OPEN) then
            call td_magnus(M, dim3, OPEN, fock_aop, F1a, F1b, rho_aop, rhonew, &
                           factorial, NBCH, dt_magnus, natom, transport_calc,  &
                           Nuc, Iz, istep, overlap, sqsm, Xmat, Ymat, Xtrans,  &
                           t/0.024190D0, M_f, fock_bop, rho_bop)
         else
            call td_magnus(M, dim3, OPEN, fock_aop, F1a, F1b, rho_aop, rhonew, &
                           factorial, NBCH, dt_magnus, natom, transport_calc,  &
                           Nuc, Iz, istep, overlap, sqsm, Xmat, Ymat, Xtrans,  &
                           t/0.024190D0, M_f)
         endif
      endif

      call g2g_timer_start('complex_rho_on_to_ao')
      call rho_aop%BChange_ONtoAO(Xtrans, M_f, .true.)
      if (OPEN) call rho_bop%BChange_ONtoAO(Xtrans, M_f, .true.)
      call g2g_timer_stop('complex_rho_on_to_ao')
      call g2g_timer_sum_pause("TD - Propagation")

      ! The real part of the density matrix in the atomic orbital basis is copied
      ! in Pmat_vec(1,2,3,...,MM) to compute the corresponding Fock matrix.
      if (OPEN) then
         call rho_aop%Gets_dataC_AO(rho_aux(:,:,1))
         call rho_bop%Gets_dataC_AO(rho_aux(:,:,2))
         call sprepack_ctr('L', M, rhoalpha, rho_aux(MTB+1:MTB+M,MTB+1:MTB+M,1))
         call sprepack_ctr('L', M, rhobeta,  rho_aux(MTB+1:MTB+M,MTB+1:MTB+M,2))

         Pmat_vec = rhoalpha + rhobeta
      else
         call rho_aop%Gets_dataC_AO(rho_aux(:,:,1))
         call sprepack_ctr('L', M, Pmat_vec, rho_aux(MTB+1:MTB+M,MTB+1:MTB+M,1))
      endif

!TBDFT: temporary TBDFT don't store restarts
      ! Stores the density matrix each 500 steps as a restart.

      call rho_aop%Gets_dataC_ON(rho(:,:,1))
      if (OPEN) call rho_bop%Gets_dataC_ON(rho(:,:,2))

      if ((writedens) .and. ( (mod(istep, td_rst_freq) == 0) .or. &
         (istep == ntdstep) )) then

         if (OPEN) then
            restart_filename='td_a.restart'
            if (istep.eq.ntdstep) restart_filename='td_last_a.restart'
            if (propagator.eq.2) then
               call write_td_restart_magnus(rho_aux(:,:,1), F1a(:,:,1), &
                                            F1b(:,:,1), M, restart_filename)
            else
               call write_td_restart_verlet(rho_aux(:,:,1), M, restart_filename)
            endif
            restart_filename='td_b.restart'
            if (istep.eq.ntdstep) restart_filename='td_last_b.restart'
            if (propagator.eq.2) then
               call write_td_restart_magnus(rho_aux(:,:,2), F1a(:,:,2), &
                                            F1b(:,:,2), M, restart_filename)
            else
               call write_td_restart_verlet(rho_aux(:,:,2), M, restart_filename)
            endif
         else
            restart_filename='td.restart'
            if (istep.eq.ntdstep) restart_filename='td_last.restart'
            if (propagator.eq.2) then
               call write_td_restart_magnus(rho_aux(:,:,1), F1a(:,:,1), &
                                            F1b(:,:,1), M, restart_filename)
            else
               call write_td_restart_verlet(rho_aux(:,:,1), M, restart_filename)
            endif
         end if
      endif

      ! Compute the trace of the density matrix for population analysis.
      if (transport_calc) call transport_rho_trace(M, dim3, rho)

      ! Dipole Moment calculation.
      call td_dipole(Pmat_vec, t, tdstep, Fx, Fy, Fz, istep, propagator, &
                     is_lpfrg, 134, 2*NCO+NUNP, r, d, Iz, pc)
      call td_population(M, natom, rho_aux(MTB+1:MTB+M,MTB+1:MTB+M,:),   &
                         Smat_initial, Nuc, Iz, OPEN, istep, propagator, &
                         is_lpfrg, td_do_pop)

      ! Population analysis.
      if (transport_calc) call transport_population(M, dim3, natom, Nuc, Iz,   &
                                                    rho_aux, overlap, sqsm,    &
                                                    propagator, is_lpfrg,      &
                                                    istep, OPEN)
      ! TD step finalization.
      call tbdft_td_output(M, dim3,rho_aux, Smat_initial, istep, Iz, natom,    &
                           Nuc, OPEN)


      if (OPEN) then
          rho_aux(:,:,1) = rho_aux(:,:,1) + rho_aux(:,:,2)
      end if

      if ( (propagator > 1) .and. (is_lpfrg) ) then
          if ( mod(istep,10) == 0 ) then
             call movieprint( natom, M, istep/10, Iz, r, dcmplx( rho_aux(:,:,1) ) )
          end if
      else
          call movieprint( natom, M, istep, Iz, r, dcmplx( rho_aux(:,:,1) ) )
      end if

      call g2g_timer_stop('TD step')
      call g2g_timer_sum_pause("TD - TD Step")

 999  continue

   ! Finalization.
   call write_energies(E1, E2, En, Ens, 0.0D0, Ex, .false., 0.0D0, 0, nsol, &
                       0.0D0, Eexact, 0.0D0, 0.0D0)
   call td_deallocate_all(F1a, F1b, fock, rho, rho_aux, rhold, rhonew, &
                          factorial, Smat_initial, Xmat, Xtrans, Ymat)
   call field_finalize()

   close(134)
   call g2g_timer_stop('TD')

   return
end subroutine TD
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Subroutines used in initialisation.

subroutine td_allocate_all(M_f,M, dim3, NBCH, propagator, F1a, F1b, fock,   &
                           rho, rho_aux, rhold, rhonew, rho_0, fock_0,sqsm, &
                           factorial, Smat_initial)
   implicit none
   integer, intent(in) :: M, NBCH, propagator, M_f, dim3
   LIODBLE, allocatable, intent(inout) :: F1a(:,:,:), F1b(:,:,:), fock(:,:,:), &
                                         sqsm(:,:), factorial(:),             &
                                         fock_0(:,:,:), Smat_initial(:,:)
   TDCOMPLEX, allocatable, intent(inout) :: rho(:,:,:), rho_aux(:,:,:),      &
                                            rhold(:,:,:), rhonew(:,:,:),     &
                                            rho_0(:,:,:)

   if ( allocated(fock)        ) deallocate(fock)
   if ( allocated(rho)         ) deallocate(rho)
   if ( allocated(rho_aux)     ) deallocate(rho_aux)
   if ( allocated(rhold)       ) deallocate(rhold)
   if ( allocated(rhonew)      ) deallocate(rhonew)
   if ( allocated(sqsm)        ) deallocate(sqsm)
   if ( allocated(factorial)   ) deallocate(factorial)
   if ( allocated(fock_0)      ) deallocate(fock_0)
   if ( allocated(rho_0)       ) deallocate(rho_0)
   if ( allocated(Smat_initial)) deallocate(Smat_initial)

   allocate(sqsm(M,M), factorial(NBCH), Smat_initial(M,M))

   allocate(fock(M_f,M_f,dim3) , rho(M_f,M_f,dim3)   , rho_aux(M_f,M_f,dim3), &
            rhold(M_f,M_f,dim3), rhonew(M_f,M_f,dim3), fock_0(M,M,dim3)     , &
            rho_0(M,M,dim3))

   if (propagator == 2) then
      if ( allocated(F1a) ) deallocate(F1a)
      if ( allocated(F1b) ) deallocate(F1b)
      allocate (F1a(M_f,M_f,dim3), F1b(M_f,M_f,dim3))
   endif

   return
end subroutine td_allocate_all

subroutine td_deallocate_all(F1a, F1b, fock, rho, rho_aux, rhold, rhonew, &
                             factorial, Smat_initial, Xmat, Xtrans, Ymat)
   use typedef_cumat, only: cumat_x, cumat_r

   implicit none
   type(cumat_r), intent(inout) :: Xmat
   type(cumat_x), intent(inout) :: Xtrans, Ymat
   LIODBLE, allocatable, intent(inout) :: F1a(:,:,:), F1b(:,:,:),   &
                                               fock(:,:,:), factorial(:),&
                                               Smat_initial(:,:)
   TDCOMPLEX   , allocatable, intent(inout) :: rho(:,:,:), rho_aux(:,:,:), &
                                               rhold(:,:,:), rhonew(:,:,:)

   call Xmat%destroy()
   call Xtrans%destroy()
   call Ymat%destroy()
   if ( allocated(fock)        ) deallocate(fock)
   if ( allocated(rho)         ) deallocate(rho)
   if ( allocated(rho_aux)     ) deallocate(rho_aux)
   if ( allocated(rhold)       ) deallocate(rhold)
   if ( allocated(rhonew)      ) deallocate(rhonew)
   if ( allocated(factorial)   ) deallocate(factorial)
   if ( allocated(F1a)         ) deallocate(F1a)
   if ( allocated(F1b)         ) deallocate(F1b)
   if ( allocated(Smat_initial)) deallocate(Smat_initial)

end subroutine td_deallocate_all

subroutine td_initialise(propagator, tdstep, NBCH, dt_lpfrg, dt_magnus, &
                         factorial)
   implicit none
   integer, intent(in)  :: propagator, NBCH
   LIODBLE , intent(in)  :: tdstep
   LIODBLE , intent(out) :: dt_lpfrg, dt_magnus, factorial(NBCH)
   integer :: icount

   select case (propagator)
   ! Initialises propagator-related parameters.
      case (1)
         dt_lpfrg = tdstep
      case (2)
         dt_magnus    = tdstep
         dt_lpfrg     = tdstep * 0.10D0
         factorial(1) = 1.0D0

         do icount = 2, NBCH
            factorial(icount) = 1.0D0 / icount
         enddo
      case default
         write(*,*) "ERROR - TD: Wrong value for propagator (td_initialise)."
         stop
   end select
end subroutine td_initialise

subroutine td_integration_setup(igrid2, igpu, atom_Z)
   implicit none
   integer, intent(in)  :: igrid2, atom_Z(:)
   integer, intent(out) :: igpu

   call g2g_timer_sum_start('TD - Exchange-correlation grid setup')
   call g2g_reload_atom_positions(igrid2, atom_Z)
   call g2g_timer_sum_stop('TD - Exchange-correlation grid setup')

   call aint_query_gpu_level(igpu)
   if (igpu.gt.1) call aint_new_step()

   return
end subroutine td_integration_setup

subroutine td_integral_1e(E1, En, E1s, Ens, MM, igpu, nsol, Pmat, Fmat, Hmat,&
                          r, pc, ntatom, natom, Smat, d, Iz)
   use faint_cpu, only: int1, intsol
   use mask_ecp , only: ECP_fock
   implicit none

   integer         , intent(in) :: MM, igpu, nsol, natom, ntatom, Iz(natom)
   LIODBLE, intent(in) :: pc(ntatom), r(ntatom,3)
   LIODBLE, intent(inout) :: Fmat(MM), Hmat(MM), E1, En, E1s, Ens
   LIODBLE, allocatable, intent(in)    :: d(:,:)
   LIODBLE, allocatable, intent(inout) :: Pmat(:), Smat(:,:)

   integer :: icount

   E1 = 0.0D0 ; En = 0.0D0
   call g2g_timer_sum_start('TD - 1-e Fock')
   call g2g_timer_sum_start('TD - Nuclear attraction')
   call int1(En, Fmat, Hmat, Smat, d, r, Iz, natom, ntatom)

   call ECP_fock(MM, Hmat)
   call g2g_timer_sum_stop('TD - Nuclear attraction')

   ! 1e terms - QMMM terms.
   if ((nsol.gt.0) .or. (igpu.ge.4)) then
      call g2g_timer_sum_start('TD - QM/MM')
      if (igpu.le.1) then
         call g2g_timer_start('intsol')
         call intsol(Pmat, Hmat, Iz, pc, r, d, natom, ntatom, E1s, Ens, .true.)
         call g2g_timer_stop('intsol')
      else
         call aint_qmmm_init(nsol, r, pc)
         call g2g_timer_start('aint_qmmm_fock')
         call aint_qmmm_fock(E1s, Ens)
         call g2g_timer_stop('aint_qmmm_fock')
      endif
         call g2g_timer_sum_stop('TD - QM/MM')
   endif

   E1=0.D0
   do icount = 1, MM
      E1 = E1 + Pmat(icount) * Hmat(icount)
   enddo
   call g2g_timer_sum_stop('TD - 1-e Fock')
   return
end subroutine td_integral_1e

subroutine td_overlap_diag(M_f, M, Smat, Xmat, Xtrans, Ymat)
   use typedef_cumat, only: cumat_r, cumat_x
   use tbdft_subs   , only: getXY_TBDFT
   use ceed_data    , only: ceed_calc, Xmat_ceed

   implicit none
   integer      , intent(in)    :: M_f, M
   LIODBLE , intent(in)    :: Smat(M,M)
   type(cumat_r), intent(inout) :: Xmat
   type(cumat_x), intent(inout) :: Xtrans, Ymat

   integer                   :: icount, jcount, LWORK, info
   LIODBLE, allocatable :: WORK(:), eigenvalues(:), X_mat(:,:),  &
                                X_trans(:,:), Y_mat(:,:), X_min(:,:), &
                                Y_min(:,:)
   TDCOMPLEX   , allocatable :: aux_mat(:,:)
   TDCOMPLEX :: liocmplx

   allocate(eigenvalues(M), X_min(M,M), Y_min(M,M), X_mat(M_f,M_f), &
            X_trans(M_f,M_f), Y_mat(M_f,M_f))

   ! Diagonalization of S matrix (LAPACK).
   X_min = Smat
   if (allocated(WORK)) deallocate(WORK)
   allocate(WORK(1))
   call dsyev('V', 'L', M, X_min, M, eigenvalues, WORK, -1, info)

   LWORK = int(WORK(1))
   deallocate(WORK)
   allocate(WORK(LWORK))
   call dsyev('V', 'L', M, X_min, M, eigenvalues, WORK, LWORK, info)

   ! Here we obtain the transformation matrices X and Y for converting from the
   ! atomic orbital basis to the molecular orbital basis (truncated during
   ! linear dependency elimination). S is the overlap matrix, s is the diagonal
   ! eigenvalue matrix of S and U is the eigenvector matrix of S:
   ! X = U s^(-1/2)
   ! Matrix X's dimension is M*3M. In the first M*M terms it contains the
   ! transformation matrices and in the other M*2M terms it contains auxiliar
   ! matrices.
   do jcount = 1, M
      if (eigenvalues(jcount) < 1.0D-6) then
         write(*,*) 'WARNING - TD: Linear dependency detected in S matrix.'
         do icount = 1, M
            X_min(icount, jcount) = 0.0D0
            Y_min(icount, jcount) = 0.0D0
         enddo
      else
         do icount = 1, M
            Y_min(icount,jcount) = X_min(icount,jcount) *sqrt(eigenvalues(jcount))
            X_min(icount,jcount) = X_min(icount,jcount) /sqrt(eigenvalues(jcount))
         enddo
      endif
   enddo

   ! CEED: X_min most be stored for CEED no matter the kind of calculation
   if (ceed_calc) call Xmat_ceed%init(M, X_min)
   ! TBDFT: Xmat and Ymat are adapted for TBDFT
   call getXY_TBDFT(M, X_min, Y_min, X_mat, Y_mat)

   ! Transposes X and stores it in Xmin temporarily.
   ! Deallocates unused arrays to avoid overloading memory.
   X_trans = transpose(X_mat)
   deallocate(eigenvalues, Y_min, X_min)

   ! Stores transformation matrices.
   call Xmat%init(M_f, X_mat)

   allocate(aux_mat(M_f,M_f))
   do icount = 1, M_f
   do jcount = 1, M_f
      aux_mat(icount,jcount) = liocmplx(X_trans(icount,jcount), 0.0D0)
   enddo
   enddo
   call Xtrans%init(M_f, aux_mat)

   do icount = 1, M_f
   do jcount = 1, M_f
      aux_mat(icount,jcount) = liocmplx(Y_mat(icount,jcount), 0.0D0)
   enddo
   enddo
   call Ymat%init(M_f, aux_mat)

   deallocate(Y_mat, aux_mat, X_trans, X_mat)

end subroutine td_overlap_diag

subroutine td_coulomb_precalc(igpu, MEMO, r, d, natom, ntatom)
   use faint_cpu, only: int3mem
   implicit none
   integer         , intent(in)    :: igpu, natom, ntatom
   LIODBLE, intent(in)    :: r(ntatom,3), d(ntatom, ntatom)
   logical         , intent(inout) :: MEMO

   if (igpu.gt.2) then
      call g2g_timer_start('Coulomb - precalc')
      call aint_coulomb_init()
      if (igpu.eq.5) MEMO = .false.
   endif
   if (MEMO) then
      call int3mem(r, d, natom, ntatom)
   endif
   call g2g_timer_stop('Coulomb - precalc')

   return
end subroutine td_coulomb_precalc

! Subroutines used in propagation.
subroutine td_get_time(t, tdstep, istep, propagator, is_lpfrg)
   implicit none
   integer, intent(in)    :: istep, propagator
   LIODBLE , intent(in)    :: tdstep
   logical, intent(in)    :: is_lpfrg
   LIODBLE , intent(inout) :: t

   select case (propagator)
      case (1)
         t = (istep-1) * tdstep
      case (2)
         if (is_lpfrg) then
            t = (istep-1) * tdstep * 0.1
         else
            t = 20 * tdstep
            t = t + (istep-200)*tdstep
         endif
      case default
   end select

   t = t * 0.024190D0
   return
end subroutine td_get_time

subroutine td_check_prop(is_lpfrg, propagator, istep, lpfrg_steps, fock_rst,&
                         verbose)
   implicit none
   integer, intent(in)  :: propagator, istep, lpfrg_steps, verbose
   logical, intent(in)  :: fock_rst
   logical, intent(out) :: is_lpfrg

   is_lpfrg = ((propagator.eq.1) .or. (((propagator.eq.2) .and. &
              (istep.lt.lpfrg_steps)) .and. (.not.fock_rst)))
   if ( (is_lpfrg).and.(istep.eq.1) ) then
      if (verbose .gt. 2) write(*,'(A)') '  TD - Starting Verlet Propagation'
   endif
   if ( (.not.(is_lpfrg)).and.(((istep-1).eq.lpfrg_steps)) ) then
      if (verbose .gt. 2) write(*,'(A)') '  TD - Starting Magnus Propagation'
   endif
   return
end subroutine td_check_prop

subroutine td_calc_energy(E, E1, E2, En, Ex, Es, Ehf, MM, Pmat, Fmat, Fmat2, &
                          Gmat, Ginv, Hmat, is_lpfrg, time, M, open_shell, r,&
                          d, natom, ntatom, MEMO, nElecs, atom_z, mm_crg)
   use faint_cpu , only: int3lu
   use field_subs, only: field_calc
   use propagators,only: do_TDexactExchange
   implicit none
   integer, intent(in)     :: MM, natom, ntatom, M, atom_z(natom), nElecs
   logical, intent(in)     :: is_lpfrg, open_shell, memo
   LIODBLE , intent(in)    :: time, r(ntatom,3), d(natom,natom), mm_crg(:)
   LIODBLE , intent(inout) :: E, E1, E2, En, Ex, Es, Ehf, Hmat(MM), Pmat(MM), &
                              Fmat(MM), Fmat2(MM), Ginv(:), Gmat(:)
   integer :: icount

   E1 = 0.0D0; E = 0.0D0
   if (is_lpfrg) then
      Ehf = 0.0d0
      call g2g_timer_sum_start("TD - Coulomb")
      call int3lu(E2, Pmat, Fmat2, Fmat, Gmat, Ginv, Hmat, open_shell,&
                  MEMO)
      call g2g_timer_sum_pause("TD - Coulomb")
      call g2g_timer_sum_start("TD - Exc")
      call g2g_solve_groups(0,Ex,0)
      call g2g_timer_sum_pause("TD - Exc")
      call g2g_timer_sum_start("TD - Exact Exchange")
      call do_TDexactExchange(Fmat,Fmat2,Ehf,MM,M,open_shell)
      call g2g_timer_sum_pause("TD - Exact Exchange")
   endif

   ! ELECTRIC FIELD CASE - Perturbation type: Gaussian (default).
   call field_calc(E1, time, Pmat, Fmat2, Fmat, r, d, natom, ntatom, open_shell, &
                   nElecs, atom_z, mm_crg)

   ! Add 1e contributions to E1.
   do icount = 1, MM
      E1 = E1 + Pmat(icount) * Hmat(icount)
   enddo
   E = E1 + E2 + En + Ex + Ehf
   if (ntatom > natom) E = E + Es

end subroutine td_calc_energy

subroutine td_dipole(rho, t, tdstep, Fx, Fy, Fz, istep, propagator, is_lpfrg, &
                     uid, nElecs, pos, dist, atom_z, mm_chrg)
   use fileio, only: write_dipole_td, write_dipole_td_header
   use properties, only: dipole
   implicit none
   integer, intent(in) :: istep, propagator, uid, nElecs, atom_z(:)
   logical, intent(in) :: is_lpfrg
   LIODBLE, intent(in) :: Fx, Fy, Fz, t, tdstep, rho(:)
   LIODBLE, intent(in) :: pos(:,:), dist(:,:), mm_chrg(:)
   LIODBLE :: dipxyz(3)

   if(istep.eq.1) then
      call write_dipole_td_header(tdstep, Fx, Fy, Fz, uid)
   endif
   if ((propagator.gt.1).and.(is_lpfrg)) then
      if (mod ((istep-1),10) == 0) then
         call g2g_timer_start('DIPOLE_TD')
         call dipole(dipxyz, rho, nElecs, pos, dist, atom_z, mm_chrg)
         call g2g_timer_stop('DIPOLE_TD')
         call write_dipole_td(dipxyz, t, uid)
      endif
   else
      call g2g_timer_start('DIPOLE_TD')
      call dipole(dipxyz, rho, nElecs, pos, dist, atom_z, mm_chrg)
      call g2g_timer_stop('DIPOLE_TD')
      call write_dipole_td(dipxyz, t, uid)
   endif

   return
end subroutine td_dipole

subroutine td_population(M, natom, rho, Smat_init, Nuc, Iz, open_shell, &
                         nstep, propagator, is_lpfrg, do_pop)
   use properties, only: print_mulliken
   use ECP_mod   , only: ecpmode, IzECP
   implicit none
   integer  , intent(in) :: M, natom, Nuc(M), Iz(natom), nstep, &
                            propagator, do_pop
   logical  , intent(in) :: open_shell, is_lpfrg
   LIODBLE  , intent(in) :: Smat_init(M,M)
   TDCOMPLEX, intent(in) :: rho(:,:,:)

   integer, allocatable :: true_iz(:)
   LIODBLE, allocatable :: real_rho(:,:), real_rho_b(:,:)

   integer          :: icount, jcount

   if (do_pop == 0) return
   if (.not. (mod(nstep, do_pop) == 0)) return
   if ((.not. (mod(nstep, do_pop*10) == 0)) .and. (propagator > 1) &
       .and. (is_lpfrg)) return
       
   allocate(true_iz(size(Iz,1)))
   true_iz = Iz
   if (ecpmode) true_iz = IzECP      

   allocate(real_rho(M,M))
   if (open_shell) then  
      allocate(real_rho_b(M,M))

      do icount = 1, M
      do jcount = 1, M
         real_rho(icount, jcount)   = real(rho(icount, jcount, 1))
         real_rho_b(icount, jcount) = real(rho(icount, jcount, 2))
      enddo
      enddo
      call print_mulliken(real_rho, real_rho_b, Smat_init, Nuc, Iz, true_iz)
      deallocate(real_rho_b)
   else
      do icount = 1, M
      do jcount = 1, M
         real_rho(icount, jcount) = real(rho(icount, jcount, 1))
      enddo
      enddo
      call print_mulliken(real_rho, Smat_init, Nuc, Iz, true_iz)
   endif
   deallocate(real_rho, true_iz)
end subroutine td_population


subroutine td_orbital_population(rho_alf, rho_bet, open_shell, nstep, &
                                 do_pop, propagator, is_lpfrg, basis_m)
   use fileio          , only: write_orbital_population
   use typedef_operator, only: operator

   implicit none
   integer       , intent(in)    :: nstep, propagator, do_pop, basis_m
   logical       , intent(in)    :: open_shell, is_lpfrg
   type(operator), intent(inout) :: rho_alf, rho_bet
  
   TDCOMPLEX , allocatable :: tmp_mat(:,:)
   LIODBLE   , allocatable :: eivec(:,:), eival(:), eivec2(:,:), eival2(:)

 
   if (do_pop == 0) return
   if (.not. (mod(nstep, do_pop) == 0)) return
   if ((.not. (mod(nstep, do_pop*10) == 0)) .and. (propagator > 1) &
        .and. (is_lpfrg)) return

   allocate(eival(basis_m), eivec(basis_m, basis_m), tmp_mat(basis_m, basis_m))

   ! Sacar dataC, poner en data, y despues diag

   call rho_alf%gets_dataC_ON(tmp_mat)
   eivec = dble(tmp_mat)
   call rho_alf%sets_data_ON(eivec)
   call rho_alf%diagon_datamat(eivec, eival)

   if (open_shell) then
      allocate(eival2(basis_m), eivec2(basis_m, basis_m))
      call rho_bet%gets_dataC_ON(tmp_mat)
      eivec2 = dble(tmp_mat)
      call rho_bet%sets_data_ON(eivec2)

      call rho_bet%diagon_datamat(eivec2, eival2)

      call write_orbital_population(eival, eival2)
   else
      call write_orbital_population(eival)
   endif

   deallocate(eival, eivec, tmp_mat)
   if (open_shell) deallocate(eival2, eivec2)
   
end subroutine td_orbital_population

subroutine td_bc_fock(M_f, M, Fmat, fock_op, Xmat, istep, time)
   use tbdft_data      , only: tbdft_calc, MTB
   use tbdft_subs      , only: chimeraTBDFT_evol
   use fockbias_subs   , only: fockbias_apply
   use typedef_operator, only: operator
   use typedef_cumat   , only: cumat_r

   implicit none
   integer       , intent(in)    :: M, M_f, istep
   LIODBLE  , intent(in)    :: time
   type(cumat_r) , intent(in)    :: Xmat
   LIODBLE  , intent(inout) :: Fmat(:)
   type(operator), intent(inout) :: fock_op

   LIODBLE, allocatable :: fock_0(:,:), fock(:,:)

   allocate(fock_0(M,M), fock(M_f,M_f))

   call g2g_timer_start('fock')
   call spunpack('L', M, Fmat, fock_0)

   if (tbdft_calc /= 0) then
      call chimeraTBDFT_evol(M, fock_0, fock, istep)
   else
      fock = fock_0
   endif
   call fockbias_apply(time, fock)

   call fock_op%Sets_data_AO(fock)
   call fock_op%BChange_AOtoON(Xmat, M_f)
   call fock_op%Gets_data_ON(fock)
   call sprepack('L', M, Fmat, fock(MTB+1:MTB+M,MTB+1:MTB+M))
   call g2g_timer_stop('fock')

   deallocate(fock_0, fock)
end subroutine td_bc_fock

subroutine td_verlet(M, M_f, dim3, OPEN, fock_aop, rhold, rho_aop, rhonew, &
                     istep, td_eu_step, Im, dt_lpfrg, transport_calc, natom,   &
                     Nuc, Iz, overlap, sqsm, Ymat, Xtrans, fock_bop, rho_bop)
   use transport_subs  , only: transport_propagate
   use tbdft_data      , only: tbdft_calc, rhold_AOTB, rhonew_AOTB, MTB
   use tbdft_subs      , only: transport_TB
   use typedef_operator, only: operator
   use typedef_cumat   , only: cumat_x
   use ceed_data       , only: ceed_calc
   use ceed_subs       , only: ceed_fock_calculation
   implicit none

   logical       , intent(in)    :: OPEN
   integer       , intent(in)    :: M, M_f, natom, Nuc(M), Iz(natom), dim3
   integer       , intent(in)    :: istep, td_eu_step
   logical       , intent(in)    :: transport_calc
   LIODBLE       , intent(in)    :: dt_lpfrg, overlap(:,:), sqsm(M,M)
   TDCOMPLEX     , intent(in)    :: Im
   type(cumat_x) , intent(in)    :: Ymat, Xtrans
   TDCOMPLEX     , intent(inout) :: rhold(M_f,M_f, dim3), rhonew(M_f,M_f, dim3)
   type(operator), intent(inout) :: fock_aop, rho_aop
   type(operator), intent(inout), optional :: fock_bop, rho_bop

   LIODBLE  , allocatable :: fock_aux(:,:,:)
   TDCOMPLEX, allocatable :: rho(:,:,:), rho_aux(:,:,:)
   TDCOMPLEX              :: liocmplx

   allocate(rho(M_f, M_f, dim3), rho_aux(M_f,M_f,dim3))
   call rho_aop%Gets_dataC_ON(rho(:,:,1))
   if (OPEN) call rho_bop%Gets_dataC_ON(rho(:,:,2))

   if (istep == 1) then
      call g2g_timer_start('cuconmut')
      call fock_aop%Commut_data_c(rho(:,:,1), rhold(:,:,1), M_f)
      if (OPEN) call fock_bop%Commut_data_c(rho(:,:,2), rhold(:,:,2), M_f)
      call g2g_timer_stop('cuconmut')
      rhold = rho + real(dt_lpfrg,COMPLEX_SIZE/2) * (Im * rhold)
   endif

   if ((transport_calc) .and. (istep >= 3)) then
      call rho_aop%Gets_dataC_AO(rho_aux(:,:,1))
      if (OPEN) call rho_bop%Gets_dataC_AO(rho_aux(:,:,2))
      call transport_propagate(M, dim3, natom, Nuc, Iz, 1, istep, overlap, &
                               sqsm, rho_aux, Ymat, OPEN)
   endif

   if (.not.OPEN) then
      call transport_TB(M, dim3, rho_aux ,Ymat, istep, OPEN, rho_aop)
   else
      call transport_TB(M, dim3, rho_aux ,Ymat, istep, OPEN, rho_aop, rho_bop)
   end if

   call g2g_timer_start('commutator')
   call fock_aop%Commut_data_c(rho(:,:,1), rhonew(:,:,1), M_f)
   if (OPEN) call fock_bop%Commut_data_c(rho(:,:,2), rhonew(:,:,2), M_f)

   !Including Euler steps if they are required
   if ((td_eu_step /= 0).and.(mod(istep, td_eu_step)==0)) then
      rhonew = rho - liocmplx(0.5d0,0.0d0)*real(dt_lpfrg,COMPLEX_SIZE/2) *     &
                     (Im * rhonew)
   else
      rhonew = rhold - real(dt_lpfrg,COMPLEX_SIZE/2) * (Im * rhonew)
   end if
   call g2g_timer_stop('commutator')

   !Including CEED term
   if(ceed_calc) then
      allocate(fock_aux(M_f,M_f,dim3))

      call fock_aop%Gets_data_ON(fock_aux(:,:,1))
      call rho_aop%Gets_dataC_ON(rho_aux(:,:,1))
      if (OPEN) then
         call fock_bop%Gets_data_ON(fock_aux(:,:,2))
         call rho_bop%Gets_dataC_ON(rho_aux(:,:,2))
      end if
      call ceed_fock_calculation(fock_aux(MTB+1:M_f,MTB+1:M_f,:),              &
                                 rho_aux(MTB+1:M_f,MTB+1:M_f,:), M, istep,     &
                                 dim3, OPEN)
      if ((td_eu_step /= 0).and.(mod(istep, td_eu_step)==0)) then
         rhonew = rhonew + real(dt_lpfrg,COMPLEX_SIZE/2) * rho_aux
      else
         rhonew = rhonew + liocmplx(2.0d0,0.0d0)*real(dt_lpfrg,COMPLEX_SIZE/2)*&
                           rho_aux
      end if

      deallocate(fock_aux)
   end if

   !Transport: Add the driving term to the propagation.
   if ((istep >= 3) .and. (transport_calc)) then
      write(*,*) 'Transport: Adding driving term to the density.'
      rhonew = rhonew - rho_aux
   endif

   ! Density update (rhold-->rho, rho-->rhonew)
   rhold = rho
   rho   = rhonew

   call rho_aop%Sets_dataC_ON(rho(:,:,1))
   if (OPEN) call rho_bop%Sets_dataC_ON(rho(:,:,2))

   ! TBDFT: rhonew and rhold in AO is store for charge calculations of TBDFT
   if (tbdft_calc == 1) then
      rhold_AOTB  = rhold
      rhonew_AOTB = rhonew
      call Xtrans%change_base(rhold_AOTB(:,:,1), 'dir')
      call Xtrans%change_base(rhonew_AOTB(:,:,1), 'dir')

      if (OPEN) then
         call Xtrans%change_base(rhold_AOTB(:,:,2), 'dir')
         call Xtrans%change_base(rhonew_AOTB(:,:,2), 'dir')
      endif
   endif

   deallocate(rho, rho_aux)
   
end subroutine td_verlet

subroutine td_magnus(M, dim3, OPEN, fock_aop, F1a, F1b, rho_aop, rhonew,       &
                     factorial, NBCH, dt_magnus, natom, transport_calc, Nuc,   &
                     Iz, istep, overlap, sqsm, Xmat, Ymat, Xtrans, time, M_f, &
                     fock_bop, rho_bop)
   use transport_subs  , only: transport_propagate
   use tbdft_data      , only: tbdft_calc, MTB, rhold_AOTB, rhonew_AOTB
   use tbdft_subs      , only: chimeraTBDFT_evol, transport_TB
   use propagators     , only: predictor, magnus
   use typedef_operator, only: operator
   use typedef_cumat   , only: cumat_r, cumat_x
   implicit none


   type(cumat_r) , intent(in) :: Xmat
   type(cumat_x) , intent(in) :: Xtrans, Ymat
   type(operator), intent(inout)           :: fock_aop, rho_aop
   type(operator), intent(inout), optional :: fock_bop, rho_bop

   logical  , intent(in)      :: transport_calc, OPEN
   integer  , intent(in)      :: dim3
   integer  , intent(in)      :: M, NBCH, istep, natom, Nuc(M), Iz(natom), M_f
   LIODBLE   , intent(in)      :: dt_magnus, factorial(NBCH), time
   LIODBLE   , intent(inout)   :: F1a(M_f,M_f, dim3), F1b(M_f,M_f, dim3),   &
                                 overlap(:,:), sqsm(M,M)
   TDCOMPLEX, intent(inout)  :: rhonew(M_f,M_f,dim3)
   TDCOMPLEX, allocatable     :: rho(:,:,:), rho_aux(:,:,:)
   LIODBLE, allocatable        :: fock_aux(:,:,:), fock(:,:,:)

   allocate(rho(M_f,M_f,dim3), rho_aux(M_f,M_f,dim3),                      &
            fock_aux(M_f,M_f, dim3), fock(M_f, M_f, dim3))

   call fock_aop%Gets_data_ON(fock(:,:,1))
   call rho_aop%Gets_dataC_ON(rho(:,:,1))

   if (OPEN) then
      call fock_bop%Gets_data_ON(fock(:,:,2))
      call rho_bop%Gets_dataC_ON(rho(:,:,2))
   end if

   if (transport_calc) then
      call rho_aop%Gets_dataC_AO(rho_aux(:,:,1))
      if (OPEN) call rho_bop%Gets_dataC_AO(rho_aux(:,:,2))
      call transport_propagate(M, dim3, natom, Nuc, Iz, 2, istep, overlap,     &
                               sqsm, rho_aux(:,:,1), Ymat, OPEN)
   endif

! TBDFT: this if is temporary, it is to conserve the atomic part of TBDFT in
! predictor.
   if (tbdft_calc /= 0) then
      call chimeraTBDFT_evol(M,fock(MTB+1:MTB+M,MTB+1:MTB+M,1), &
                             fock_aux(:,:,1), istep)

      if (OPEN) call chimeraTBDFT_evol(M,fock(MTB+1:MTB+M,MTB+1:MTB+M,2), &
                                       fock_aux(:,:,2), istep)

      ! TBDFT: rhold in AO is store for charge calculations of TBDFT
      if (tbdft_calc == 1) then
         rhold_AOTB = rho
         call Xtrans%change_base(rhold_AOTB(:,:,1), 'dir')
         if (OPEN) call Xtrans%change_base(rhold_AOTB(:,:,2), 'dir')
      end if

      fock = fock_aux
   endif

!DLVN_TBDFT: calculating the drving term
   if (.not.OPEN) then
      call transport_TB(M, dim3, rho_aux, Ymat, istep, OPEN, rho_aop)
   else
      call transport_TB(M, dim3, rho_aux, Ymat, istep, OPEN, rho_aop, rho_bop)
   end if


   call g2g_timer_start('predictor')
   call predictor(F1a, F1b, fock, rho, factorial, Xmat, Xtrans, dt_magnus, &
                  time, M_f, MTB, dim3)
   call g2g_timer_stop('predictor')
   call g2g_timer_start('magnus')
   call magnus(fock(:,:,1), rho(:,:,1), rhonew(:,:,1), M_f, NBCH, dt_magnus,  &
               factorial)
   if (OPEN) call magnus(fock(:,:,2), rho(:,:,2), rhonew(:,:,2), M_f, NBCH,   &
                         dt_magnus, factorial)
   call g2g_timer_stop('magnus')

   ! Transport: Add the driving term to the propagation.
   if (transport_calc) then
      write(*,*) 'Transport: Adding driving term to the density.'
      rhonew = rhonew - rho_aux
   endif

   ! DLVN_TBDFT: Adding driving term to the propagation
   if (tbdft_calc == 3) then
      write(*,*) 'Transport TB: Adding driving term to the density.'
      rhonew = rhonew - rho_aux
   endif

   ! TBDFT: rhonew in AO is store for charge calculations of TBDFT
   if (tbdft_calc == 1) then
      rhonew_AOTB = rhonew
      call Xtrans%change_base(rhonew_AOTB(:,:,1), 'dir')
      if (OPEN) call Xtrans%change_base(rhonew_AOTB(:,:,2), 'dir')
   endif

   ! Density update and Fock storage.
   F1a = F1b
   F1b = fock
   rho = rhonew
   call fock_aop%Sets_data_ON(fock(:,:,1))
   call rho_aop%Sets_dataC_ON(rho(:,:,1))

   if (OPEN) then
      call fock_bop%Sets_data_ON(fock(:,:,2))
      call rho_bop%Sets_dataC_ON(rho(:,:,2))
   end if

   return
end subroutine td_magnus

subroutine calc_trace(matrix, msize, message)

   implicit none
   integer         , intent(in) :: msize
   character(len=*), intent(in) :: message
   LIODBLE, intent(in) :: matrix(msize, msize)

   integer          :: icount
   LIODBLE :: trace

   trace = 0.0D0

   do icount = 1, msize
      trace = trace + matrix(icount, icount)
   enddo

   write(*,*) "Trace of ", message, " equals to ", trace
   return
end subroutine calc_trace

subroutine calc_trace_c(matrix, msize, message)

   implicit none
   integer         , intent(in) :: msize
   character(len=*), intent(in) :: message
   integer          :: icount
   TDCOMPLEX, intent(in) :: matrix(msize, msize)
   TDCOMPLEX :: trace
   TDCOMPLEX :: liocmplx
   trace = liocmplx(0.0D0, 0.0D0)

   do icount = 1, msize
      trace = trace + matrix(icount, icount)
   enddo

   write(*,*) "Trace of ", message, " equals to ", trace
   return
end subroutine calc_trace_c
end module time_dependent
