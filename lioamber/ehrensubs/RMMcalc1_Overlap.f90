!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine RMMcalc1_Overlap(Smat,Energy)
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use maskrmm     , only: rmmget_fock
  use garcha_mod  , only: M
  use faint_cpu77 , only: int1
  implicit none
  real*8,intent(out) :: Smat(M,M)
  real*8,intent(out) :: Energy
  integer            :: igpu

  call g2g_timer_start('RMMcalc1')
  call aint_query_gpu_level(igpu)
  if (igpu.gt.1) call aint_new_step()
  call int1(Energy)
  call rmmget_fock(Smat)
  call g2g_timer_stop('RMMcalc1')


end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
