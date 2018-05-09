!_o%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrendyn_main( energy_o, dipmom_o )
!------------------------------------------------------------------------------!
!
!  stored_densM1 and stored_densM2 are stored in ON basis, except for the
!  first step
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use garcha_mod, &
   &  only: M, natom, atom_mass, nucpos, nucvel, qm_forces_ds, qm_forces_total &
   &      , first_step, propagator

   use td_data, &
   &  only: tdstep

   use ehrendata, &
   &  only: stored_time, stored_energy, stored_dipmom                          &
   &      , stored_densM1, stored_densM2                                       &
   &      , rsti_funit, rsto_funit, nustep_count, elstep_count                 &
   &       , ndyn_steps, edyn_steps, wdip_nfreq, wdip_fname                    &
   &      , rsti_loads, rsti_fname, rsto_saves, rsto_nfreq, rsto_fname

   implicit none
   real*8,intent(inout) :: dipmom_o(3), energy_o
   real*8               :: dipmom(3)  , energy  , energy0

   real*8  :: time, dtn, dte, dtaux
   integer :: elstep_local, elstep_keeps
   integer :: substep, substeps
   integer :: nn, kk

   logical :: first_nustep
   logical :: load_restart
   logical :: rhomid_in_ao
   logical :: missing_last

   real*8, allocatable, dimension(:,:) :: nucfor_ds
   real*8, allocatable, dimension(:,:) :: Smat, Sinv
   real*8, allocatable, dimension(:,:) :: Lmat, Umat, Linv, Uinv
   real*8, allocatable, dimension(:,:) :: Fock, Fock0
   real*8, allocatable, dimension(:,:) :: Bmat, Dmat

   complex*16, allocatable, dimension(:,:) :: RhoOld, RhoMid, RhoNew
   complex*16, allocatable, dimension(:,:) :: RhoMidF
   complex*16, allocatable, dimension(:,:) :: Tmat

   logical, parameter :: velocity_recalc = .true.
!
!
!
!  Preliminaries
!------------------------------------------------------------------------------!
   call g2g_timer_start('ehrendyn - nuclear step')
!  BEWARE OF COUNTING => EXTRA STEP WHEN NOT DOING RESTART...
!   if (first_step) return
   nustep_count = nustep_count + 1
   time = stored_time

   allocate( nucfor_ds(3,natom) )
   allocate( Smat(M,M), Sinv(M,M) )
   allocate( Lmat(M,M), Umat(M,M), Linv(M,M), Uinv(M,M) )
   allocate( Fock(M,M), Fock0(M,M) )
   allocate( RhoOld(M,M), RhoMid(M,M), RhoNew(M,M), RhoMidF(M,M) )
   allocate( Bmat(M,M), Dmat(M,M), Tmat(M,M) )

   dtn = tdstep
   dte = ( tdstep / edyn_steps )

   first_nustep = (nustep_count == 1)
   load_restart = (first_nustep).and.(rsti_loads)
   rhomid_in_ao = (first_nustep).and.(.not.rsti_loads)
   missing_last = (first_nustep).and.(.not.rsti_loads)

   if (first_nustep) stored_energy = energy_o
   if (load_restart) call ehrenaux_rsti( rsti_fname, &
   &  natom, qm_forces_total, nucvel, M, stored_densM1, stored_densM2 )
!
!
!  Update velocities, calculate fixed fock, load last step dens matrices
!------------------------------------------------------------------------------!
   if (velocity_recalc) then
      dtaux = dtn/2.0d0 - dte/2.0d0
      call ehrenaux_updatevel( natom, atom_mass, qm_forces_total, nucvel, dtaux )
   else
      call ehrenaux_updatevel( natom, atom_mass, qm_forces_total, nucvel, dtn )
   endif

   energy0 = 0.0d0
   call RMMcalc0_Init()
   call RMMcalc1_Overlap( Smat, energy0 )
   call ehrenaux_cholesky( M, Smat, Lmat, Umat, Linv, Uinv, Sinv )
   call RMMcalc2_FockMao( Fock0, energy0 )

   RhoOld = stored_densM1
   RhoMid = stored_densM2
   if (rhomid_in_ao) then
      RhoMid = matmul(RhoMid, Lmat)
      RhoMid = matmul(Umat, RhoMid)
      stored_densM2 = RhoMid
   endif
!
!
!
!  ELECTRONIC STEP CYCLE
!------------------------------------------------------------------------------!
   elstep_keeps = ceiling( real(edyn_steps) / 2.0 )

   do elstep_local = 1, edyn_steps
      call g2g_timer_start('ehrendyn - electronic step')
      elstep_count = elstep_count + 1
      dipmom(:) = 0.0d0
      energy = energy0
      Fock = Fock0

      if (velocity_recalc) call ehrenaux_updatevel &
      &  ( natom, atom_mass, qm_forces_total, nucvel, dte )

      call ehrendyn_step( missing_last, propagator, time, dte, M, natom,       &
                        & nucpos, nucvel, nucfor_ds, Sinv, Uinv, Linv,         &
                        & RhoOld, RhoMid, RhoNew, Fock, dipmom, energy )

      RhoOld = RhoMid
      RhoMid = RhoNew

      if ( elstep_local == elstep_keeps ) qm_forces_ds = nucfor_ds
      time = time + dte * 0.0241888d0

      call ehrenaux_writedip(elstep_count, wdip_nfreq, stored_time, dipmom,    &
      &    wdip_fname)

      if (rsto_saves) call ehrenaux_rsto( &
      &  rsto_fname, rsto_nfreq, ndyn_steps*edyn_steps, elstep_count,          &
      & natom, qm_forces_total, nucvel, M, RhoOld, RhoMid )

      call g2g_timer_stop('ehrendyn - electronic step')

   enddo
!
!
!
!  Finalizations
!------------------------------------------------------------------------------!
   stored_densM1 = RhoOld
   stored_densM2 = RhoMid

   dipmom_o = stored_dipmom
   energy_o = stored_energy
   stored_dipmom = dipmom
   stored_energy = energy
   stored_time = time

   deallocate( Smat, Sinv )
   deallocate( Lmat, Umat, Linv, Uinv )
   deallocate( Fock, Fock0 )
   deallocate( RhoOld, RhoMid, RhoNew, RhoMidF )
   deallocate( Bmat, Dmat, Tmat )
   call g2g_timer_stop('ehrendyn - nuclear step')

end subroutine ehrendyn_main
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
