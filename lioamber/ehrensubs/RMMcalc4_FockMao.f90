!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine RMMcalc4_FockMao( DensMao, FockMao, DipMom, Energy )
!------------------------------------------------------------------------------!
!
! DESCRIPTION (Time is in ps?fs?)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use maskrmm

   use faint_cpu  , only: int2, int3mem, int3lu, intsol, intfld

   use field_data , only: epsilon, a0

   use garcha_mod, &
   &only: M, Md, RMM, kkind, kkinds, cool, cools, igrid2                       &
       &, natom, Iz, NCO, Nunp, total_time, d, ntatom, r, open, pc

   use ehrendata, &
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
   integer  :: MM, MMd, igpu, M7, M9, M3, M5, M11
   logical  :: MEMO

!  For electric field application
   real*8   :: FieldNow(3)
   real*8   :: factor, g, Qc
   real*8   :: dip_times_field, strange_term
   real*8   :: field_shape, time_fact, time_dist, laser_freq

   MM=M*(M+1)/2
   MMd=Md*(Md+1)/2
   M3=1+MM ! Pew
   M5=M3+MM ! now S, also F later
   M7  = 1 + 3*MM
   M9  = M7 + MMd
   M11=M9+MMd ! Hmat
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

   idx0=3*MM+2*MMd
   if (igpu.le.1) then
      call intsol(RMM(1:MM), RMM(idx0+1:idx0+MM+1), Iz, pc, r, d, natom, &
                  ntatom, Energy_SolvF, Energy_SolvT, .true.)
   else
      call aint_qmmm_fock(Energy_SolvF,Energy_SolvT)
   endif

   call int2(RMM(M7:M7+MMd), RMM(M9:M9+MMd), r, d, ntatom)
   if (igpu.gt.2) call aint_coulomb_init()
   MEMO = .true.
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
   call int3lu(Energy_Coulomb, RMM(1:MM), RMM(M3:M3+MM), RMM(M5:M5+MM),        &
               RMM(M7:M7+MMd), RMM(M9:M9+MMd), RMM(M11:M11+MMd), open)
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
     call intfld(RMM(M3:M3+MM), RMM(M5:M5+MM), r, d, Iz, natom, ntatom, open, &
                 g, FieldNow(1), FieldNow(2), FieldNow(3))

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
