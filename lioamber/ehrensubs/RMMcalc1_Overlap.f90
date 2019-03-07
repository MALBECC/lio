!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine RMMcalc1_Overlap(Ovlap,Energy)
!
! DESCRIPTION
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use garcha_mod  , only: RMM, Smat, d, r, Iz, natom, ntatom, Fmat_vec
  use basis_data  , only: M, Md
  use faint_cpu   , only: int1
  implicit none
  real*8,intent(out) :: Ovlap(M,M)
  real*8,intent(out) :: Energy
  integer            :: igpu,  MM, M11
  MM = M*(M+1)/2 ; M11 = 1 + MM*2+MM+Md*(Md+1)


  call g2g_timer_start('RMMcalc1')
  call aint_query_gpu_level(igpu)
  if (igpu.gt.1) call aint_new_step()
  call int1(Energy, Fmat_vec, RMM(M11:M11+MM), Smat, d, r, Iz, natom, &
            ntatom)
  call rmmget_fock(Ovlap)
  call g2g_timer_stop('RMMcalc1')


end subroutine RMMcalc1_Overlap
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
