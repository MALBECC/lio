!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Calculates the forces in the QM region.
subroutine dft_get_qm_forces(dxyzqm)
   use garcha_mod , only: natom, ntatom, nsol, r, d, Iz, first_step,   &
                          cubegen_only, number_restr, doing_ehrenfest, &
                          qm_forces_ds, qm_forces_total, Pmat_en_wgt,  &
                          Pmat_vec
   use ehrendata  , only: nullify_forces
   use faint_cpu  , only: int1G, intSG, int3G
   use fileio_data, only: verbose
   use fileio     , only: write_force_log
   use dftd3      , only: dftd3_gradients
   implicit none
   double precision, intent(out) :: dxyzqm(3,natom)
   double precision, allocatable :: ff1G(:,:),ffSG(:,:),ff3G(:,:), ffvdw(:,:)
   integer            :: igpu, katm, icrd
   double precision   :: f_r ! For restraints

   if (cubegen_only) return
   call g2g_timer_sum_start('Forces')
   allocate(ff1G(natom,3), ffSG(natom,3), ff3G(natom,3), ffvdw(natom,3))
   ff1G = 0.0D0 ; ffSG = 0.0D0 ; ff3G = 0.0D0; ffvdw = 0.0D0

   ! 1e gradients.
   call g2g_timer_start('int1G')
   call aint_query_gpu_level(igpu)
   if (igpu .lt. 4) then
      call g2g_timer_sum_start('Nuclear attraction gradients')
      call int1G(ff1G, Pmat_vec, d, r, Iz, natom, ntatom)
      call g2g_timer_sum_stop('Nuclear attraction gradients')
   elseif (nsol .le. 0) then
      call g2g_timer_sum_start('Nuclear attraction gradients')
      call int1G(ff1G, Pmat_vec, d, r, Iz, natom, ntatom)
      call aint_qmmm_forces(ff1G,0)
      call g2g_timer_sum_stop('Nuclear attraction gradients')
   endif
   call g2g_timer_stop('int1G')

   ! Overlap gradients.
   call g2g_timer_start('intSG')
   call g2g_timer_sum_start('Overlap gradients')
   if (doing_ehrenfest) then
      ffSG = -transpose(qm_forces_ds)
   else
      call intSG(ffSG, Pmat_en_wgt, r, d, natom, ntatom)
   endif
   call g2g_timer_sum_stop('Overlap gradients')
   call g2g_timer_stop('intSG')

   ! 2e and 3e gradients.
   call g2g_timer_start('int3G')
   call g2g_timer_sum_start('Coulomb+Exchange-correlation')
   call int3G(ff3G, .true., Pmat_vec, r, d, natom, ntatom)
   call g2g_timer_stop('int3G')
   call g2g_timer_sum_stop('Coulomb+Exchange-correlation')

   ! DFTD3 gradients
   call g2g_timer_sum_start("DFTD3 Gradients")
   call dftd3_gradients(ffvdw, r, natom)
   call g2g_timer_sum_stop("DFTD3 Gradients")

   ! Gets total
   do katm = 1, natom
   do icrd = 1, 3
      dxyzqm(icrd,katm) = ff1G(katm,icrd) + ffSG(katm,icrd) + ff3G(katm,icrd) + &
                          ffvdw(katm,icrd)
   enddo
   enddo

   ! FFR - Ehrenfest needs to keep track of forces
   if ( nullify_forces ) dxyzqm(:,:) = 0.0D0
   if ( doing_ehrenfest ) then
      qm_forces_total = qm_forces_ds
      qm_forces_total = qm_forces_total - transpose(ff1G)
      qm_forces_total = qm_forces_total - transpose(ff3G)
   endif

   ! Distance restrains
   if (number_restr .gt. 0) call get_restrain_forces(dxyzqm, f_r)

   ! FFR: force calculation should be separated from force passing and
   !      force writing. All can be in the same module, but different
   !      subroutines.
   if (verbose .gt. 4) call write_force_log(dxyzqm, ff1G, ffSG, ff3G, natom, &
                                            3242, first_step)

   if (nsol.le.0) call g2g_timer_sum_stop('Forces')

   ! FFR: No other place for this to go right now.
   if ( first_step ) first_step = .false.

   deallocate(ff1G, ffSG, ff3G, ffvdw)
end subroutine dft_get_qm_forces
