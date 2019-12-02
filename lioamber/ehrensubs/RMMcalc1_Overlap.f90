!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine RMMcalc1_Overlap(Ovlap,Energy)
!
! DESCRIPTION
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use garcha_mod  , only: Smat, d, r, Iz, natom, ntatom, Fmat_vec, Hmat_vec
  use basis_data  , only: M
  use faint_cpu   , only: int1
  implicit none
  real*8,intent(out) :: Ovlap(M,M)
  real*8,intent(out) :: Energy
  integer            :: igpu

  call g2g_timer_start('RMMcalc1')
  call aint_query_gpu_level(igpu)
  if (igpu.gt.1) call aint_new_step()
  call int1(Energy, Fmat_vec, Hmat_vec, Smat, d, r, Iz, natom, ntatom)
  call rmmget_fock(Ovlap)
  call g2g_timer_stop('RMMcalc1')


end subroutine RMMcalc1_Overlap
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
