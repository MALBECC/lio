!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmcalc2_focknuc( fock_mao, energy_1e, energy_solvT )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use faint_cpu  , only: int2, int3mem, intsol
   use garcha_mod , only: M, Md, RMM, kkind, kkinds, cool, cools, igrid2, MEMO,&
                          d, r, ntatom, Iz, pc, natom
   use ECP_mod    , only: ecpmode, term1e, VAAA, VAAB, VBAC, &
                        & FOCK_ECP_read, FOCK_ECP_write

   implicit none
   real*8, intent(out)    :: fock_mao(M,M)
   real*8, intent(out)    :: energy_1e
   real*8, intent(out)    :: energy_solvT

   real*8  :: energy_solvF
   integer :: kk, idx0
   integer :: MM, MMd, igpu, M9, M7
!
!
!  Initializations
!------------------------------------------------------------------------------!
   call g2g_timer_start('rmmcalc2-init')
   call rmmput_fock( fock_mao )

!  Pointer M11
   MM   = M  * (M+1)  / 2
   MMd  = Md * (Md+1) / 2
   idx0 = 3*MM + 2*MMd
   M7  = 1 + 3*MM
   M9  = M7 + MMd

   if (allocated(kkind))  deallocate(kkind)
   if (allocated(kkinds)) deallocate(kkinds)
   if (allocated(cool))   deallocate(cool)
   if (allocated(cools))  deallocate(cools)

   call g2g_reload_atom_positions(igrid2)
   call aint_query_gpu_level(igpu)
   if (igpu.gt.1) call aint_new_step()
   call g2g_timer_stop('rmmcalc2-init')
!
!
! Calculate fixed-parts of fock
!------------------------------------------------------------------------------!
   call g2g_timer_start('rmmcalc2-sol2coul')
   if (igpu.le.1) then
      call intsol(RMM(1:MM), RMM(idx0:idx0+MM), Iz, pc, r, d, natom, ntatom, &
                  energy_solvF, energy_solvT, .true.)
   else
      call aint_qmmm_fock( energy_solvF, energy_solvT )
   endif

   call int2(RMM(M7:M7+MMd), RMM(M9:M9+MMd), r, d, ntatom)
   if (igpu.gt.2) call aint_coulomb_init()
   if (igpu.eq.5) MEMO = .false.
   call g2g_timer_stop('rmmcalc2-sol2coul')

   if (MEMO) then
      call g2g_timer_start('rmmcalc2-int3mem')
      call int3mem()
      call g2g_timer_stop('rmmcalc2-int3mem')
   endif
!
!
! Calculate Pseudo-Potentials
!------------------------------------------------------------------------------!
   if (ecpmode) then
      call g2g_timer_start('rmmcalc2-ECP')

      if (FOCK_ECP_read) then
!        Variable allocation and data read from ECP_restart
         call intECP(0)

      else
!        Variable allocation and calculation of the N center terms by different
!        calls to subroutine intECP(N).
         call intECP(1)
         call intECP(2)
         call intECP(3)

      end if

      if (FOCK_ECP_write) call WRITE_ECP()
      call WRITE_POST(1)

      write(*,*) "Modifying Fock Matrix with ECP terms"
      do kk = 1, MM
!        Backups 1e terms and modifies fock
         term1e( kk )   = RMM( idx0+kk )
         RMM( idx0+kk ) = RMM( idx0+kk ) + VAAA(kk) + VAAB(kk) + VBAC(kk)
      enddo

      call g2g_timer_stop('rmmcalc2-ECP')
   end if
!
!
!  Prepare Outputs
!------------------------------------------------------------------------------!
   call g2g_timer_start('rmmcalc2-exit')

   energy_1e = 0.0d0
   do kk = 1, MM
      energy_1e = energy_1e + RMM(kk)*RMM(idx0+kk)
   enddo

   call rmmget_fock( fock_mao )
   call g2g_timer_stop('rmmcalc2-exit')

end subroutine rmmcalc2_focknuc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
