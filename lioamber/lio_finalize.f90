subroutine lio_finalize()
! Deallocation and finalizations.
   use garcha_mod , only: dipole, Smat, RealRho, sqsm, Eorbs, Eorbs_b, &
                          MO_coef_at, MO_coef_at_b, r, v, rqm, Em, Rm, &
                          pc, Iz, d, Fmat_vec, Fmat_vec2, Pmat_vec,    &
                          Hmat_vec, Ginv_vec, Gmat_vec, Pmat_en_wgt
   use ECP_mod       , only: ecpmode
   use fileio        , only: io_finish_outputs
   use basis_subs    , only: basis_deinit
   use converger_subs, only: converger_finalise
   use dftd3         , only: dftd3_finalise
 
   implicit none
   call basis_deinit() ! Deallocates basis variables.
   call io_finish_outputs(dipole, 69)
   if (ecpmode) call generalECP(4) ! Deallocates ECP variables.

   ! Deallocates global variables.
   if (allocated(Smat))         deallocate(Smat)
   if (allocated(RealRho))      deallocate(RealRho)
   if (allocated(sqsm))         deallocate(sqsm)
   if (allocated(Eorbs))        deallocate(Eorbs)
   if (allocated(Eorbs_b))      deallocate(Eorbs_b)
   if (allocated(MO_coef_at))   deallocate(MO_coef_at)
   if (allocated(MO_coef_at_b)) deallocate(MO_coef_at_b)


   deallocate(Fmat_vec, Fmat_vec2, Pmat_vec, Hmat_vec, Ginv_vec, &
              Gmat_vec, Pmat_en_wgt)

   deallocate(r, v, rqm, Em, Rm, pc, Iz, d)

   ! GPU code finalization.
   call aint_deinit()
   call g2g_timer_summary()
   call g2g_deinit()
   call converger_finalise()
   call dftd3_finalise()
end subroutine
