!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine ehrendyn(Energy,DipMom)
!------------------------------------------------------------------------------!
!
! RhoSaveA and RhoSaveB are stored in ON basis, except for the first step
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use garcha_mod, only: M, natom, tdstep, total_time, first_step, atom_mass    &
                      , nucpos, nucvel, qm_forces_ds, qm_forces_total
  implicit none
  real*8,intent(inout) :: Energy,DipMom(3)

  real*8,allocatable,dimension(:,:)     :: Smat,Sinv,Lmat,Umat,Linv,Uinv
  real*8,allocatable,dimension(:,:)     :: Fock,FockInt
  complex*16,allocatable,dimension(:,:) :: RhoOld,RhoMid,RhoNew

  real*8,allocatable,dimension(:,:)     :: Bmat,Dmat
  complex*16,allocatable,dimension(:,:) :: Tmat
  real*8                                :: dt
  real*8                                :: ux,uy,uz
  integer                               :: ii,jj,kk,idx

! Preliminaries
!------------------------------------------------------------------------------!
  print*,'Doing ehrenfest!'
  step_number = step_number + 1
  call g2g_timer_start('ehrendyn step')
  allocate( Smat(M,M), Sinv(M,M) )
  allocate( Lmat(M,M), Umat(M,M), Linv(M,M), Uinv(M,M) )
  allocate( Fock(M,M), FockInt(M,M) )
  allocate( RhoOld(M,M), RhoMid(M,M), RhoNew(M,M) )
  allocate( Bmat(M,M), Dmat(M,M), Tmat(M,M) )

  if (.not.allocated(qm_forces_total)) then
    allocate(qm_forces_total(3,natom))
    qm_forces_total=0.0d0
  endif

  if (.not.allocated(qm_forces_ds)) then
    allocate(qm_forces_ds(3,natom))
    qm_forces_ds=0.0d0
  endif

  dt=tdstep

  if ( first_step .and. restart_dyn ) then
     open( unit=rstinp_unit, file="ehren_rstinp" )
     print*,'Using restart'
     call rstload( rstinp_unit, Natom, qm_forces_total, M, RhoSaveA, RhoSaveB )
     close( unit=rstinp_unit )
  end if


! Update velocities
!------------------------------------------------------------------------------!
  do ii=1,natom
  do kk=1,3
    nucvel(kk,ii)=nucvel(kk,ii)+(1.5d0)*dt*qm_forces_total(kk,ii)/atom_mass(ii)
  enddo
  enddo

! Nuclear Force Calculation (works in AO)
!------------------------------------------------------------------------------!
  Energy=0.0d0
  call RMMcalc0_Init()
  call RMMcalc1_Overlap(Smat,Energy)
  call ehren_cholesky(M,Smat,Lmat,Umat,Linv,Uinv,Sinv)

! Esto deja la Rho correcta en RMM, pero habria que ordenarlo mejor
  RhoMid=RhoSaveB
  if (.not.first_step) then
    RhoMid=matmul(RhoMid,Linv)
    RhoMid=matmul(Uinv,RhoMid)
  endif
  call RMMcalc2_FockMao(RhoMid,Fock,DipMom,Energy)
  call calc_forceDS(natom,M,nucpos,nucvel,RhoMid,Fock,Sinv,Bmat,qm_forces_ds)



! Density Propagation (works in ON)
!------------------------------------------------------------------------------!
  Fock=matmul(Fock,Uinv)
  Fock=matmul(Linv,Fock)
  Dmat=calc_Dmat(M,Linv,Uinv,Bmat)
  Tmat=DCMPLX(Fock)+DCMPLX(0.0d0,1.0d0)*DCMPLX(Dmat)

  RhoOld=RhoSaveA
  RhoMid=RhoSaveB
  if (first_step) then
    RhoMid=matmul(RhoMid,Lmat)
    RhoMid=matmul(Umat,RhoMid)
    call ehren_verlet_e(M,-(dt/2.0d0),Tmat,RhoMid,RhoMid,RhoOld)
  endif
  call ehren_verlet_e(M,dt,Tmat,RhoOld,RhoMid,RhoNew)
  RhoSaveA=RhoMid
  RhoSaveB=RhoNew


! Saving restart
!------------------------------------------------------------------------------!
  if ( modulo(step_number,save_lapse) == 1 ) then
     save_step = .true.
  else
     save_step = .false.
  end if

  if ( step_number == last_step ) save_step = .true.

  if ( save_step ) then
     open( unit=rstout_unit, file="ehren_rstout" )
     call rstsave( rstout_unit, Natom, qm_forces_total, M, RhoSaveA, RhoSaveB )
     close( unit=rstout_unit )
  end if

! Calculation of the dipole moment
!------------------------------------------------------------------------------!
  if (first_step) then
    open(unit=134,file='x.dip')
    open(unit=135,file='y.dip')
    open(unit=136,file='z.dip')
    write(134,*) '#Time (fs) vs DIPOLE MOMENT, X COMPONENT (DEBYES)'
    write(135,*) '#Time (fs) vs DIPOLE MOMENT, Y COMPONENT (DEBYES)'
    write(136,*) '#Time (fs) vs DIPOLE MOMENT, Z COMPONENT (DEBYES)'
    total_time=0.0d0
  endif

  if (.not.first_step) then
   call dip(ux,uy,uz)
   write(134,901) total_time,ux
   write(135,901) total_time,uy
   write(136,901) total_time,uz
   print*,''
   print*,' Timer: ',total_time
   print*,''
   total_time=total_time+dt*0.0241888d0
  endif

  call g2g_timer_stop('ehrendyn step')

 901 format(F15.9,2x,F15.9)

  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
