!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine RMMcalc2_FockMao( FockMao, Energy )
!------------------------------------------------------------------------------!
!
! DESCRIPTION (Time is in ps?fs?)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use faint_cpu  , only: int2, int3mem, intsol

   use garcha_mod,  only: RMM, igrid2, MEMO, r, d, ntatom, Iz, pc, natom, &
                          Gmat_vec, Ginv_vec, Hmat_vec
   use basis_data,  only: M, Md, kkind, kkinds, cool, cools, nshelld, MM, MMd

   implicit none
   real*8,intent(out)    :: FockMao(M,M)
   real*8,intent(out)    :: Energy

   real*8   :: Energy_1e
   real*8   :: Energy_Coulomb
   real*8   :: Energy_SolvT,Energy_SolvF

   integer  :: kk, igpu, M7
!
!
!  Initializations
!------------------------------------------------------------------------------!
   call g2g_timer_start('RMMcalc2-init')
   call rmmput_fock(FockMao)

   if (allocated(kkind))  deallocate(kkind)
   if (allocated(kkinds)) deallocate(kkinds)
   if (allocated(cool))   deallocate(cool)
   if (allocated(cools))  deallocate(cools)

   call g2g_reload_atom_positions(igrid2)
   call aint_query_gpu_level(igpu)
   if (igpu.gt.1) call aint_new_step()
   call g2g_timer_stop('RMMcalc2-init')
   M7  = 1 + 3*MM
!
!
! Calculate fixed-parts of fock
!------------------------------------------------------------------------------!
   call g2g_timer_start('RMMcalc2-sol2coul')
   if (igpu.le.1) then
      call intsol(RMM(1:MM), Hmat_vec, Iz, pc, r, d, natom, &
                  ntatom, Energy_SolvF, Energy_SolvT, .true.)
   else
      call aint_qmmm_fock(Energy_SolvF,Energy_SolvT)
   endif

   call int2(RMM(M7:M7+MMd), Gmat_vec, r, d, ntatom)
   if (igpu.gt.2) call aint_coulomb_init()
   if (igpu.eq.5) MEMO = .false.
   call g2g_timer_stop('RMMcalc2-sol2coul')

   if (MEMO) then
      call g2g_timer_start('RMMcalc2-int3mem')
      call int3mem(r, d, natom, ntatom)
      call g2g_timer_stop('RMMcalc2-int3mem')
   endif
!
!
!  Prepare Outputs
!------------------------------------------------------------------------------!
   call g2g_timer_start('RMMcalc2-exit')
   Energy_1e=0.0d0
   do kk=1,MM
      Energy_1e = Energy_1e + RMM(kk) * Hmat_vec(kk)
   enddo

!  Energy=0.0d0
   Energy=Energy+Energy_1e
   Energy=Energy+Energy_Coulomb
   Energy=Energy+Energy_SolvT

   call rmmget_fock(FockMao)
   call g2g_timer_stop('RMMcalc2-exit')

end subroutine RMMcalc2_FockMao
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
