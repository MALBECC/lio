!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine RMMcalc3_FockMao( DensMao, ElecField, FockMao, DipMom, Energy )
!------------------------------------------------------------------------------!
!
! DESCRIPTION (Time is in ps?fs?)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use maskrmm,     only: rmmput_dens, rmmput_fock, rmmget_fock

   use faint_cpu  , only: intfld, int3lu

   use garcha_mod,  only: M, natom, Iz, NCO, Nunp, total_time, RMM, Md, open, &
                          r, d, ntatom

   use field_data,  only: epsilon, a0

   use ehrendata, only: eefld_on

   implicit none
   complex*16,intent(in) :: DensMao(M,M)
   real*8,intent(in)     :: ElecField(3)
   real*8,intent(inout)  :: FockMao(M,M)
   real*8,intent(inout)  :: DipMom(3)
   real*8,intent(inout)  :: Energy

   real*8   :: Energy_Coulomb
   real*8   :: Energy_Exchange
   real*8   :: Energy_Efield
   integer  :: kk, MM, MMd, M3, M5, M7, M9, M11

!  For electric field
   real*8   :: factor, g, Qc
   real*8   :: dip_times_field, strange_term
!
!
!  Calculate unfixed Fock in RMM - int3lu and solve_groups
!------------------------------------------------------------------------------!
   MM =M *(M+1)/2
   MMd=Md*(Md+1)/2
   M3=1+MM ! Pew
   M5=M3+MM ! now S, also F later
   M7=M5+MM ! G matrix
   M9=M7+MMd ! G inverted
   M11=M9+MMd ! Hmat

   call g2g_timer_start('RMMcalc3-solve3lu')
   call rmmput_fock( FockMao )
   call rmmput_dens( DensMao )
   call int3lu(Energy_Coulomb, RMM(1:MM), RMM(M3:M3+MM), RMM(M5:M5+MM),        &
               RMM(M7:M7+MMd), RMM(M9:M9+MMd), RMM(M11:M11+MMd), open)
   call g2g_solve_groups( 0, Energy_Exchange, 0 )
   call g2g_timer_stop('RMMcalc3-solve3lu')
!
!
!  Calculate unfixed Fock in RMM - electric field
!------------------------------------------------------------------------------!
   call g2g_timer_start('RMMcalc3-field')
   call dip( DipMom(1), DipMom(2), DipMom(3) )
   write(666,*) eefld_on
   if (eefld_on) then
      g = 1.0d0
      factor = 2.54d0
      Qc = (-2.0d0) * NCO+Nunp
      do kk = 1, natom
         Qc = Qc + Iz(kk)
      end do

      call intfld(RMM(M3:M3+MM), RMM(M5:M5+MM), r, d, Iz, natom, ntatom, open, &
                  g, ElecField(1), ElecField(2), ElecField(3) )

      dip_times_field = 0.0d0
      dip_times_field = dip_times_field + ElecField(1) * DipMom(1)
      dip_times_field = dip_times_field + ElecField(2) * DipMom(2)
      dip_times_field = dip_times_field + ElecField(3) * DipMom(3)
      strange_term = (0.5d0) * (1.0d0 - 1.0d0/epsilon) * Qc**2 / a0

      Energy_Efield = 0.0d0
      Energy_Efield = Energy_Efield - g * dip_times_field / factor
      Energy_Efield = Energy_Efield - strange_term

      write(666,*)
      write(666,*) DipMom(1), DipMom(2), DipMom(3)
      write(666,*) ElecField(1), ElecField(2), ElecField(3)
      write(666,*) Energy_Efield

   end if
   call g2g_timer_stop('RMMcalc3-field')
!
!
!  Prepare outputs
!------------------------------------------------------------------------------!
   call g2g_timer_start('RMMcalc3-exit')
!  Energy = 0.0d0
   Energy = Energy + Energy_Coulomb
   Energy = Energy + Energy_Exchange
   Energy = Energy + Energy_Efield

   call rmmget_fock( FockMao )
   call g2g_timer_stop('RMMcalc3-exit')

end subroutine RMMcalc3_FockMao
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
