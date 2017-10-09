!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine RMMcalc1_Overlap(Smat,Energy)
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use maskrmm
  use garcha_mod, only:M
  implicit none
  real*8,intent(out) :: Smat(M,M)
  real*8,intent(out) :: Energy
  integer            :: igpu

  call aint_query_gpu_level(igpu)
  if (igpu.gt.1) call aint_new_step()
  call int1(Energy)
  call rmmget_fock(Smat)

  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
