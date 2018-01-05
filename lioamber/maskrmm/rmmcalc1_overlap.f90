!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmcalc1_overlap( Smat, energy_nuc )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use faint_cpu77 , only: int1
   use garcha_mod  , only: M

   implicit none
   real*8 , intent(out) :: Smat(M,M)
   real*8 , intent(out) :: energy_nuc
   integer              :: igpu

   call g2g_timer_start('rmmcalc1')
   call aint_query_gpu_level( igpu )
   if (igpu.gt.1) call aint_new_step()
   call int1( energy_nuc )
   call rmmget_fock( Smat )
   call g2g_timer_stop('rmmcalc1')

end subroutine rmmcalc1_overlap
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
