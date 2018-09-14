!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine RMMcalc1_Overlap(Ovlap,Energy)
!
! DESCRIPTION
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use maskrmm     , only: rmmget_fock
  use garcha_mod  , only: RMM, Smat, Nuc, a, c, d, r, Iz, ncont, NORM, natom, M, Md, nshell, ntatom
  use faint_cpu   , only: int1
  implicit none
  real*8,intent(out) :: Ovlap(M,M)
  real*8,intent(out) :: Energy
  integer            :: igpu,  MM, M5, M11
  MM = M*(M+1)/2 ; M5 = 1 + MM*2; M11 = M5+MM+Md*(Md+1)


  call g2g_timer_start('RMMcalc1')
  call aint_query_gpu_level(igpu)
  if (igpu.gt.1) call aint_new_step()
  call int1(Energy,RMM(M5:M5+MM), RMM(M11:M11+MM), Smat, d, r, Iz, natom, &
            ntatom)
  call rmmget_fock(Ovlap)
  call g2g_timer_stop('RMMcalc1')


end subroutine RMMcalc1_Overlap
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
