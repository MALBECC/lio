!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmcalc3_fockele( dens_mao, elec_field, uses_field,                 &
                           & fock_mao, dipole, energy_coul, energy_xc,         &
                           & energy_field, energy_ecp )
!
!  Time is in ps?fs?
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use garcha_mod  , only: M, RMM, natom, Iz, NCO, Nunp, Md, open, r, d, ntatom
   use field_data  , only: a0, epsilon
   use faint_cpu   , only: int3lu, intfld
   use ECP_mod     , only: ecpmode, VAAA, VAAB, VBAC

   implicit none
   complex*16, intent(in)    :: dens_mao(M,M)
   real*8    , intent(in)    :: elec_field(3)
   logical   , intent(in)    :: uses_field
   real*8    , intent(inout) :: fock_mao(M,M)
   real*8    , intent(inout) :: dipole(3)
   real*8    , intent(inout) :: energy_coul
   real*8    , intent(inout) :: energy_xc
   real*8    , intent(inout) :: energy_field
   real*8    , intent(inout) :: energy_ecp

   integer :: MM, kk, MMd, M3, M5, M7, M9, M11
   real*8  :: factor, g, Qc
   real*8  :: dip_times_field, strange_term
!
!
!  Calculate unfixed Fock in RMM - int3lu and solve_groups
!------------------------------------------------------------------------------!
   MM=M*(M+1)/2
   MMd=Md*(Md+1)/2
   M3=1+MM ! Pew
   M5=M3+MM ! now S, also F later
   M7=M5+MM ! G matrix
   M9=M7+MMd ! G inverted
   M11=M9+MMd ! Hmat

   call g2g_timer_start('rmmcalc3-solve3lu')
   call rmmput_fock( fock_mao)
   call rmmput_dens( dens_mao )
   call int3lu( energy_coul, RMM(1:MM), RMM(M3:M3+MM), RMM(M5:M5+MM), &
                RMM(M7:M7+MMd), RMM(M9:M9+MMd), RMM(M11:M11+MMd), open)
   call g2g_solve_groups( 0, energy_xc, 0 )
   call g2g_timer_stop('rmmcalc3-solve3lu')
!
!
!  Calculate unfixed Fock in RMM - electric field
!------------------------------------------------------------------------------!
   call g2g_timer_start('rmmcalc3-dipole')
   call dip( dipole(1), dipole(2), dipole(3) )
   call g2g_timer_stop('rmmcalc3-dipole')

   if (uses_field) then
      call g2g_timer_start('rmmcalc3-field')
      g = 1.0d0
      factor = 2.54d0

      Qc = (-2.0d0) * NCO+Nunp
      do kk = 1, natom
         Qc = Qc + Iz(kk)
      end do

      call intfld(RMM(M3:M3+MM), RMM(M5:M5+MM), r, d, Iz, natom, ntatom, open, &
                  g, elec_field(1), elec_field(2), elec_field(3))

      dip_times_field = 0.0d0
      dip_times_field = dip_times_field + elec_field(1) * dipole(1)
      dip_times_field = dip_times_field + elec_field(2) * dipole(2)
      dip_times_field = dip_times_field + elec_field(3) * dipole(3)
      strange_term = (0.5d0) * (1.0d0 - 1.0d0/epsilon) * Qc**2 / a0

      energy_field = 0.0d0
      energy_field = energy_field - g * dip_times_field / factor
      energy_field = energy_field - strange_term

      call g2g_timer_stop('rmmcalc3-field')
   endif
!
!
!  Prepare outputs
!------------------------------------------------------------------------------!
   call g2g_timer_start('rmmcalc3-exit')

   energy_ecp = 0.0d0
   if (ecpmode) then
      MM = M * (M+1) / 2
      do kk = 1, MM
         energy_ecp = energy_ecp + RMM(kk)*(VAAA(kk)+VAAB(kk)+VBAC(kk))
      enddo
   end if

   call rmmget_fock( fock_mao )
   call g2g_timer_stop('rmmcalc3-exit')

end subroutine rmmcalc3_fockele
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
