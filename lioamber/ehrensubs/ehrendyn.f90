!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrendyn( energy_o, dipmom_o )
!------------------------------------------------------------------------------!
!
!  RhoSaveA and RhoSaveB are stored in ON basis, except for the first step
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use garcha_mod, &
   &  only: M, natom, tdstep, total_time, first_step, atom_mass                &
   &      , nucpos, nucvel, qm_forces_ds, qm_forces_total

   use lionml_data, &
   &  only: ndyn_steps, edyn_steps &
   &      , rsti_loads, rsti_fname, rsto_saves, rsto_nfreq, rsto_fname

   use ehrendata, &
   &  only: RhoSaveA, RhoSaveB, rsti_funit, rsto_funit, step_number            &
   &      , StoredEnergy

   implicit none
   real*8,intent(inout) :: energy_o, dipmom_o(3)
   real*8               :: energy,   dipmom(3)

   real*8  :: dipmom_norm
   real*8  :: nucvel_update(3)
   real*8  :: dtn
   integer :: nn, kk
   logical :: first_from_scratch

   real*8,allocatable,dimension(:,:)     :: Smat, Sinv
   real*8,allocatable,dimension(:,:)     :: Lmat, Umat, Linv, Uinv
   real*8,allocatable,dimension(:,:)     :: Fock, Fock0
   complex*16,allocatable,dimension(:,:) :: RhoOld, RhoMid, RhoNew
   real*8,allocatable,dimension(:,:)     :: Bmat, Dmat
   complex*16,allocatable,dimension(:,:) :: Tmat
!
!
!
!  Preliminaries
!------------------------------------------------------------------------------!
   call g2g_timer_start('ehrendyn step')
   print*,'Doing ehrenfest!'
   step_number = step_number + 1
   allocate( Smat(M,M), Sinv(M,M) )
   allocate( Lmat(M,M), Umat(M,M), Linv(M,M), Uinv(M,M) )
   allocate( Fock(M,M), Fock0(M,M) )
   allocate( RhoOld(M,M), RhoMid(M,M), RhoNew(M,M) )
   allocate( Bmat(M,M), Dmat(M,M), Tmat(M,M) )

   dtn=tdstep

   if (first_step) then
   if (rsti_loads) then
      call ehrenrsti_load( rsti_fname, rsti_funit, Natom, qm_forces_total,  &
                         & nucvel, M, RhoSaveA, RhoSaveB )
   endif
   endif

   first_from_scratch = (first_step).and.(.not.rsti_loads)
!
!
!
!  Update velocities
!------------------------------------------------------------------------------!
   do nn=1,natom
   do kk=1,3
      nucvel_update(kk) = (1.5d0) * dtn * qm_forces_total(kk,nn) / atom_mass(nn)
      nucvel(kk,nn)     = nucvel(kk,nn) + nucvel_update(kk)
   enddo
   enddo
!
!
!
!  Nuclear Force Calculation (works in AO)
!------------------------------------------------------------------------------!
   energy=0.0d0
   call RMMcalc0_Init()
   call RMMcalc1_Overlap(Smat,energy)
   call ehren_cholesky(M,Smat,Lmat,Umat,Linv,Uinv,Sinv)

!  Esto deja la Rho correcta en RMM, pero habria que ordenarlo mejor
   RhoMid=RhoSaveB
   if (.not.first_from_scratch) then
      RhoMid=matmul(RhoMid,Linv)
      RhoMid=matmul(Uinv,RhoMid)
   endif
   call RMMcalc2_FockMao(Fock,energy)
   call RMMcalc3_FockMao(RhoMid,Fock,dipmom,energy)
   call calc_forceDS(natom,M,nucpos,nucvel,RhoMid,Fock,Sinv,Bmat,qm_forces_ds)
!
!
!
!  Density Propagation (works in ON)
!------------------------------------------------------------------------------
   call g2g_timer_start('ehrendyn - density propagation')
   Fock=matmul(Fock,Uinv)
   Fock=matmul(Linv,Fock)
   Dmat=calc_Dmat(M,Linv,Uinv,Bmat)
   Tmat=DCMPLX(Fock)+DCMPLX(0.0d0,1.0d0)*DCMPLX(Dmat)

   RhoOld=RhoSaveA
   RhoMid=RhoSaveB
   if (first_from_scratch) then
      RhoMid=matmul(RhoMid,Lmat)
      RhoMid=matmul(Umat,RhoMid)
      call ehren_verlet_e(M,-(dtn/2.0d0),Tmat,RhoMid,RhoMid,RhoOld)
   endif
   call ehren_verlet_e(M,dtn,Tmat,RhoOld,RhoMid,RhoNew)
   RhoSaveA=RhoMid
   RhoSaveB=RhoNew
   call g2g_timer_stop('ehrendyn - density propagation')
!
!
!
! Calculation of the dipole moment (TODO: REMOVE?)
!------------------------------------------------------------------------------!
   if (first_step) then
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
         & step_number, Natom, qm_forces_total, nucvel, M, RhoSaveA, RhoSaveB)
   endif

   dipmom_o = dipmom
   energy_o = StoredEnergy
   StoredEnergy = energy

   deallocate( Smat, Sinv )
   deallocate( Lmat, Umat, Linv, Uinv )
   deallocate( Fock, Fock0 )
   deallocate( RhoOld, RhoMid, RhoNew )
   deallocate( Bmat, Dmat, Tmat )
   call g2g_timer_stop('ehrendyn step')

901 format(F15.9,2x,F15.9)
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
