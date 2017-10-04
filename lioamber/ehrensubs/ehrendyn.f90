!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrendyn( energy_o, dipmom_o )
!------------------------------------------------------------------------------!
!
!  RhoSaveA and RhoSaveB are stored in ON basis, except for the first step
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use garcha_mod, &
   &  only: M, natom, tdstep, total_time, atom_mass, nucpos, nucvel            &
   &      , qm_forces_ds, qm_forces_total

   use lionml_data, &
   &  only: ndyn_steps, edyn_steps, propagator &
   &      , rsti_loads, rsti_fname, rsto_saves, rsto_nfreq, rsto_fname

   use ehrendata, &
   &  only: StoredEnergy, RhoSaveA, RhoSaveB, rsti_funit, rsto_funit           &
   &      , nustep_count, elstep_count

   implicit none
   real*8,intent(inout) :: dipmom_o(3), energy_o
   real*8               :: dipmom(3)  , energy  , energy0
   real*8               :: dipmom_norm

   real*8  :: nucvel_update(3), time_factor
   real*8  :: dtn, dte, dtaux
   integer :: elstep_local, elstep_keeps
   integer :: substep, substeps
   integer :: nn, kk

   logical :: first_nustep
   logical :: load_restart
   logical :: rhomid_in_ao
   logical :: missing_last

   real*8, allocatable, dimension(:,:) :: kept_forces
   real*8, allocatable, dimension(:,:) :: Smat, Sinv
   real*8, allocatable, dimension(:,:) :: Lmat, Umat, Linv, Uinv
   real*8, allocatable, dimension(:,:) :: Fock, Fock0
   real*8, allocatable, dimension(:,:) :: Bmat, Dmat

   complex*16, allocatable, dimension(:,:) :: RhoOld, RhoMid, RhoNew
   complex*16, allocatable, dimension(:,:) :: RhoMidF
   complex*16, allocatable, dimension(:,:) :: Tmat
!
!
!
!  Preliminaries
!------------------------------------------------------------------------------!
   print*,'Doing ehrenfest!'

   call g2g_timer_start('ehrendyn - nuclear step')
   nustep_count = nustep_count + 1

   allocate( kept_forces(3,natom) )
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

   if (load_restart) then
      call ehrenrsti_load( rsti_fname, rsti_funit, natom, qm_forces_total,  &
                         & nucvel, M, RhoSaveA, RhoSaveB )
   endif

!
!
!
!  Update velocities, calculate fixed fock, load last step dens matrices
!------------------------------------------------------------------------------!
   do nn=1,natom
   do kk=1,3
!      time_factor       = (1.0d0) * dtn - (0.5d0) * dte
!      nucvel_update(kk) = time_factor * qm_forces_total(kk,nn) / atom_mass(nn)
      nucvel_update(kk) = dtn * qm_forces_total(kk,nn) / atom_mass(nn)
      nucvel(kk,nn)     = nucvel(kk,nn) + nucvel_update(kk)
   enddo
   enddo

   energy0 = 0.0d0
   call RMMcalc0_Init()
   call RMMcalc1_Overlap( Smat, energy0 )
   call ehren_cholesky( M, Smat, Lmat, Umat, Linv, Uinv, Sinv )
   call RMMcalc2_FockMao( Fock0, energy0 )

   RhoOld = RhoSaveA
   RhoMid = RhoSaveB
   if (rhomid_in_ao) then
      RhoMid   = matmul(RhoMid, Lmat)
      RhoMid   = matmul(Umat, RhoMid)
      RhoSaveB = RhoMid
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
      substeps = 20

      if (missing_last) then
         dtaux = (-dte) / ( (2.0d0)*(substeps) )
         RhoOld = RhoMid
         call ehrenstep( 1, dtaux, M, natom, nucpos, nucvel, qm_forces_ds,     &
                       & Sinv, Uinv, Linv, RhoOld, RhoMid, RhoNew, Fock,       &
                       & energy, dipmom )
         RhoOld = RhoNew
         dtaux = (dte) / (substeps)

         do substep = 1, substeps
            dipmom(:) = 0.0d0
            energy = energy0
            Fock = Fock0
            call ehrenstep( 1, dtaux, M, natom, nucpos, nucvel, qm_forces_ds,  &
                          & Sinv, Uinv, Linv, RhoOld, RhoMid, RhoNew, Fock,    &
                          & energy, dipmom )
            RhoOld = RhoMid
            RhoMid = RhoNew
         enddo
         RhoOld = RhoSaveB
         missing_last = .false.

      else
         call ehrenstep( propagator, dte, M, natom, nucpos, nucvel,            &
                       & qm_forces_ds, Sinv, Uinv, Linv, RhoOld, RhoMid,       &
                       & RhoNew, Fock, energy, dipmom )
         RhoOld = RhoMid
         RhoMid = RhoNew

      end if

      do nn = 1, natom
      do kk = 1, 3
         nucvel_update(kk) = dte * qm_forces_total(kk,nn) / atom_mass(nn)
         nucvel(kk,nn)     = nucvel(kk,nn) + nucvel_update(kk)
      enddo
      enddo

      if ( elstep_local == elstep_keeps ) kept_forces = qm_forces_ds
      call g2g_timer_stop('ehrendyn - electronic step')

   enddo

   RhoSaveA = RhoOld
   RhoSaveB = RhoMid
   qm_forces_ds = kept_forces
!
!
!
! Calculation of the dipole moment (TODO: REMOVE?)
!------------------------------------------------------------------------------!
   if (first_nustep) then
      call write_dipole(dipmom, 0, 134, .true.)
      total_time = 0.0d0
   else
      dipmom_norm = sqrt( dipmom(1)**2 + dipmom(2)**2 + dipmom(3)**2 )
      call write_dipole( dipmom, dipmom_norm, 134, .false.)  
      total_time = total_time + dtn * 0.0241888d0
   endif
!
!
!
!  Finalizations
!------------------------------------------------------------------------------!
   if (rsto_saves) then
      call ehrenrsto_save( rsto_fname, rsto_funit, rsto_nfreq, ndyn_steps,     &
         & nustep_count, Natom, qm_forces_total, nucvel, M, RhoSaveA, RhoSaveB)
   endif

   dipmom_o = dipmom
   energy_o = StoredEnergy
   StoredEnergy = energy

   deallocate( Smat, Sinv )
   deallocate( Lmat, Umat, Linv, Uinv )
   deallocate( Fock, Fock0 )
   deallocate( RhoOld, RhoMid, RhoNew, RhoMidF )
   deallocate( Bmat, Dmat, Tmat )
   call g2g_timer_stop('ehrendyn - nuclear step')

901 format(F15.9,2x,F15.9)
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
