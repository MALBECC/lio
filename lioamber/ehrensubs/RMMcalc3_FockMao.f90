!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine RMMcalc3_FockMao( DensMao, ElecField, FockMao, DipMom, Energy )
!------------------------------------------------------------------------------!
!
! DESCRIPTION (Time is in ps?fs?)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use faint_cpu  , only: intfld, int3lu

   use garcha_mod,  only: natom, Iz, NCO, Nunp, Pmat_vec, open, &
                          r, d, ntatom, MEMO, Fmat_vec, Fmat_vec2, Ginv_vec, &
                          Hmat_vec, Gmat_vec, pc
   use properties, only: dipole
   
   use basis_data, only: M

   use field_data,  only: epsilon, a0

   use ehrendata, only: eefld_on

   implicit none
   complex(kind=8),intent(in) :: DensMao(M,M)
   LIODBLE,intent(in)     :: ElecField(3)
   LIODBLE,intent(inout)  :: FockMao(M,M)
   LIODBLE,intent(inout)  :: DipMom(3)
   LIODBLE,intent(inout)  :: Energy

   LIODBLE   :: Energy_Coulomb
   LIODBLE   :: Energy_Exchange
   LIODBLE   :: Energy_Efield
   integer  :: kk

!  For electric field
   LIODBLE   :: factor, g, Qc
   LIODBLE   :: dip_times_field, strange_term
!
!
!  Calculate unfixed Fock in RMM - int3lu and solve_groups
   call g2g_timer_start('RMMcalc3-solve3lu')
   call rmmput_fock( FockMao )
   call rmmput_dens( DensMao )
   call int3lu(Energy_Coulomb, Pmat_vec, Fmat_vec2, Fmat_vec, Gmat_vec, &
               Ginv_vec, Hmat_vec, open, MEMO)
   call g2g_solve_groups( 0, Energy_Exchange, 0 )
   call g2g_timer_stop('RMMcalc3-solve3lu')
!
!
!  Calculate unfixed Fock in RMM - electric field
!------------------------------------------------------------------------------!
   call g2g_timer_start('RMMcalc3-field')
   call dipole( DipMom, Pmat_vec, 2*NCO+Nunp, r, d, Iz, pc, .true.)
   if (eefld_on) then
      write(666,*) eefld_on
      g = 1.0d0
      factor = 2.54d0
      Qc = (-2.0d0) * NCO+Nunp
      do kk = 1, natom
         Qc = Qc + Iz(kk)
      end do

      call intfld(Fmat_vec2, Fmat_vec, r, d, natom, ntatom, open, &
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
