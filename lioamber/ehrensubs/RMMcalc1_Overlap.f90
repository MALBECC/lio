!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine RMMcalc1_Overlap(Ovlap,Energy)
!
! DESCRIPTION
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use maskrmm     , only: rmmget_fock
  use garcha_mod  , only: RMM, Smat, Nuc, a, c, d, r, Iz, ncont, NORM, natom, M, Md
  use faint_cpu   , only: int1
  implicit none
  real*8,intent(out) :: Ovlap(M,M)
  real*8,intent(out) :: Energy
  integer            :: igpu

  call g2g_timer_start('RMMcalc1')
  call aint_query_gpu_level(igpu)
  if (igpu.gt.1) call aint_new_step()
!  call int1(Energy)
  call int1(Energy,RMM,Smat,Nuc,a,c,d,r,Iz,ncont,NORM,natom,M,Md)
  call rmmget_fock(Ovlap)
  call g2g_timer_stop('RMMcalc1')


end subroutine RMMcalc1_Overlap
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
