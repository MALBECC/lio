!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine RMMcalc4_FockMao( DensMao, FockMao, DipMom, Energy )
!------------------------------------------------------------------------------!
!
! DESCRIPTION (Time is in ps?fs?)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use maskrmm
   use faint_cpu77, only: intsol, int2, int3mem, int3lu, intfld

   use field_data , only: epsilon, a0

   use garcha_mod, &
   &only: M, Md, RMM, kkind, kkinds, cool, cools, igrid2                       &
       &, natom, Iz, NCO, Nunp, total_time

   use lionml_data, &
   &only: eefld_on, eefld_ampx, eefld_ampy, eefld_ampz, eefld_wavelen          &
       &, eefld_timegih, eefld_timegfh, eefld_timepos, eefld_timeamp

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

   integer  :: kk, idx0
   integer  :: MM, MMd, igpu
   logical  :: MEMO

!  For electric field application
   real*8   :: FieldNow(3)
   real*8   :: factor, g, Qc
   real*8   :: dip_times_field, strange_term
   real*8   :: field_shape, time_fact, time_dist, laser_freq
!
!
! Calculate fixed-parts of fock
!------------------------------------------------------------------------------!
   call g2g_timer_start('RMMcalc4')
   call g2g_timer_start('RMMcalc4-start')
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
   call g2g_timer_stop('RMMcalc4-start')
   if (MEMO) then
      call g2g_timer_start('RMMcalc4-int3mem')
      call int3mem()
      call g2g_timer_stop('RMMcalc4-int3mem')
   endif

!
!  Calculate unfixed Fock in RMM - int3lu and solve_groups
!------------------------------------------------------------------------------!
   call g2g_timer_start('RMMcalc4-solve3lu')
   call rmmput_dens(DensMao)
   call int3lu(Energy_Coulomb)
   call g2g_solve_groups(0,Energy_Exchange,0)
   call g2g_timer_stop('RMMcalc4-solve3lu')

!
!  Calculate unfixed Fock in RMM - electric field
!------------------------------------------------------------------------------!
   call g2g_timer_start('RMMcalc4-field')
   if (eefld_on) then
      g = 1.0d0
      factor = 2.54d0
      field_shape = 1.0d0

      Qc = (-2.0d0) * NCO+Nunp
      do kk = 1, natom
         Qc = Qc + Iz(kk)
      end do

!     Time Gaussian Shape
      if ( ( (eefld_timegih).and.( total_time < eefld_timepos ) ) &
      &.or.( (eefld_timegfh).and.( total_time > eefld_timepos ) ) ) then
         time_fact = (-1.0d0) / (eefld_timeamp)
         time_dist = (total_time - eefld_timepos)**2
         field_shape = field_shape * exp( (time_fact) * (time_dist) )
      end if

!     Laser shape
      if (eefld_wavelen > 0.0d0) then
         laser_freq = (6.28318530718d0) * (299.792458d0) / (eefld_wavelen)
!        [freq fs-1]=       [2pi]       *  [c in nm/fs]  /   [long in nm]
         field_shape = field_shape * sin( laser_freq * total_time )
      end if

!    Apply field
     FieldNow(1) = eefld_ampx * field_shape
     FieldNow(2) = eefld_ampy * field_shape
     FieldNow(3) = eefld_ampz * field_shape
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

   end if
   call g2g_timer_stop('RMMcalc4-field')
!
!
! Calculate Energy
!------------------------------------------------------------------------------!
   call g2g_timer_start('RMMcalc4-exit')
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
!------------------------------------------------------------------------------!
   call rmmget_fock(FockMao)
   call g2g_timer_stop('RMMcalc4-exit')
   call g2g_timer_stop('RMMcalc4')

end subroutine RMMcalc4_FockMao
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
