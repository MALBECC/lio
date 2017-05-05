!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine RMMcalc2_FockMao( DensMao, FockMao, DipMom, Energy )
!
! Time is in ps
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use maskrmm
  use garcha_mod, only:M,Md,RMM,kkind,kkinds,cool,cools,igrid2, &
                     & total_time,field,Fx,Fy,Fz,epsilon,a0,    &
                     & natom,tdstep,Iz,NCO,Nunp
  implicit none
  complex*16,intent(in) :: DensMao(M,M)
  real*8,intent(out)    :: FockMao(M,M)
  real*8,intent(out)    :: DipMom(3)
  real*8,intent(out)    :: Energy

  real*8   :: Energy_1e
  real*8   :: Energy_Coulomb
  real*8   :: Energy_SolvT,Energy_SolvF
  real*8   :: Energy_Exchange
  real*8   :: Energy_Efield

  real*8   :: FieldNow(3)
  real*8   :: factor, g, Qc
  real*8   :: dip_times_field, strange_term
  integer  :: kk,idx0
  integer  :: MM,MMd,igpu
  logical  :: MEMO

  logical  :: external_field_is_present
  logical  :: laser_is_on
  real*8   :: laser_freq, laser_long
  real*8   :: alpha, timeof_pert, timeof_peak
  real*8   :: efield_shape, shapeshift
!
!
! Calculate fixed-parts of fock
!--------------------------------------------------------------------!
  if (allocated(kkind))  deallocate(kkind)
  if (allocated(kkinds)) deallocate(kkinds)
  if (allocated(cool))   deallocate(cool)
  if (allocated(cools))  deallocate(cools)

  call g2g_reload_atom_positions(igrid2)
  call aint_query_gpu_level(igpu)
  if (igpu.gt.1) call aint_new_step()

  if (igpu.le.1) then
    call intsol(Energy_SolvF,Energy_SolvT,.true.)
  else
    call aint_qmmm_fock(Energy_SolvF,Energy_SolvT)
  endif

  call int2()

  if (igpu.gt.2) call aint_coulomb_init()
  if (igpu.eq.5) MEMO = .false.
  if (MEMO) then
    call g2g_timer_start('int3mem')
    call int3mem()
    call g2g_timer_stop('int3mem')
  endif

!
! Calculate unfixed Fock in RMM - int3lu and solve_groups
!--------------------------------------------------------------------!
  call rmmput_dens(DensMao)
  call g2g_timer_start('g2g-solve + int3lu')
  call int3lu(Energy_Coulomb)
  call g2g_solve_groups(0,Energy_Exchange,0)
  call g2g_timer_stop('g2g-solve + int3lu')

!
! Calculate unfixed Fock in RMM - electric field
!--------------------------------------------------------------------!
!  if ( field ) external_field_is_present = .true.
  external_field_is_present = .true.

  if ( external_field_is_present ) then
    g = 1.0d0
    factor = 2.54d0
    efield_shape = 1.0d0

    Qc=-2*NCO+Nunp
    do kk=1,natom
      Qc=Qc+Iz(kk)
    end do

!!! Time Gaussian Shape (time in fs)
!    timeof_peak =  50.0d0 * ( tdstep * 0.0241888d0 )
!    timeof_pert = 100.0d0 * ( tdstep * 0.0241888d0 )
!    alpha       =   0.2d0 / ( tdstep * 0.0241888d0 )**2

    timeof_peak = 20.0d0
    timeof_pert = 50.0d0
!    alpha       = 0.1d-7 / ( tdstep * 0.0241888d0 )**2
    alpha       = 0.1d-7 / ( 2.06706875d-002 * 0.0241888d0 )**2

    shapeshift   = exp( (-alpha) * ( total_time - timeof_peak )**2 )
    efield_shape = efield_shape * shapeshift
!    if ( total_time <= timeof_pert ) then
!       efield_shape = efield_shape * shapeshift
!    end if

!!! Laser shape
    laser_is_on = .true.
!    laser_long = 139.3d0 !nm
    laser_long = 83.22d0 !nm
!    laser_long =123.22d0 !nm
    laser_freq = (6.28318530718d0) * (299.792458d0) / (laser_long)
!   [freq fs-1] = [2pi] * [c in nm/fs] / [long]
    shapeshift = sin( laser_freq * total_time )
    if ( laser_is_on ) then
       efield_shape = efield_shape * shapeshift
    endif

    FieldNow(1) = Fx * efield_shape
    FieldNow(2) = Fy * efield_shape
    FieldNow(3) = Fz * efield_shape

    call dip( DipMom(1), DipMom(2), DipMom(3) )
    call intfld( g, FieldNow(1), FieldNow(2), FieldNow(3) )
    dip_times_field = 0.0d0
    dip_times_field = dip_times_field + FieldNow(1) * DipMom(1)
    dip_times_field = dip_times_field + FieldNow(2) * DipMom(2)
    dip_times_field = dip_times_field + FieldNow(3) * DipMom(3)
    strange_term = (0.5d0) * (1.0d0 - 1.0d0/epsilon) * Qc**2 / a0

    Energy_Efield = 0.0d0
    Energy_Efield = Energy_Efield - g * dip_times_field / factor
    Energy_Efield = Energy_Efield - strange_term

    write(666,*) total_time, efield_shape, Energy_Efield
  end if


!
! Calculate Energy
!--------------------------------------------------------------------!
  MM=M*(M+1)/2
  MMd=Md*(Md+1)/2
  idx0=3*MM+2*MMd
  Energy_1e=0.0d0
  do kk=1,MM
    Energy_1e=Energy_1e+RMM(kk)*RMM(idx0+kk)
  enddo

!  Energy=0.0d0
  Energy=Energy+Energy_1e
  Energy=Energy+Energy_Coulomb
  Energy=Energy+Energy_SolvT
  Energy=Energy+Energy_Exchange
  Energy=Energy+Energy_Efield

!
! Extract FockMao from RMM
!--------------------------------------------------------------------!
  call rmmget_fock(FockMao)


  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
