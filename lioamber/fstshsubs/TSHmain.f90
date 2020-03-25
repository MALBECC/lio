subroutine TSHmain(CoefA, EneA, Etot)
! This routine perform TSH dynamic between differents excited
! states. 
! This routine do not perform electronic interpolation
! The electronic interpolation is performed after the 
! velocities actualization. This is called by HYBRID 
use garcha_mod  , only: OPEN, NCO, PBE0, r
use excited_data, only: lresp, nstates, libint_recalc, fittExcited
use excitedsubs , only: fcaApp, basis_initLR, fca_restored, linear_response, &
                        basis_deinitLR
use basis_data  , only: M, c_raw, c, a
use fstsh_data  , only: call_number, C_scf_old, WFcis_old, all_states, &
                        tsh_nucStep, sigma_old, Sovl_old, Sovl_now, tsh_file, &
                        a_old, c_old, r_old
   implicit none

   LIODBLE, intent(in) :: CoefA(:,:), EneA(:)
   LIODBLE, intent(inout) :: Etot

   integer :: NCOlr, Mlr, Nvirt, Ndim, ndets
   LIODBLE, allocatable :: Xexc(:,:), Eexc(:)
   LIODBLE, allocatable :: WFcis(:,:), sigma(:,:)
   LIODBLE, allocatable :: C_scf(:,:), E_scf(:)

   if (OPEN  .eqv. .true. ) then
      print*, "TSH doesn't work in Open shell yet"
      stop
   endif

   ! Initialization of Libint
   if ( (.not. PBE0) .and.  (.not. fittExcited) ) then
      call g2g_timer_sum_start('Libint init')
      call g2g_libint_init(c_raw,libint_recalc)
      call g2g_timer_sum_stop('Libint init')
   endif

   ! Linear Response Calculation
   !   We calculate all excited states involved in the dynamic
   call fcaApp(CoefA,EneA,C_scf,E_scf,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
   call basis_initLR(C_scf,M,Mlr,NCOlr,Nvirt)
   call g2g_saverho( )
   allocate(Xexc(Ndim,nstates),Eexc(nstates))
   call linear_response(C_scf,E_scf,Xexc,Eexc,M,Mlr,Nvirt,NCOlr,Ndim,0)
   call fca_restored(CoefA,EneA,C_scf,E_scf,Xexc,M,Mlr,Nvirt,NCO,NCOlr,&
                     Ndim,nstates)

   ! Set Cis Wavefunctions dimensions
   ndets = (NCO * Nvirt)*2 + 1 ! ES + GS determinants
   if ( tsh_nucStep == 0 ) then
      allocate(WFcis_old(ndets,all_states),sigma_old(all_states,all_states))
   endif

   if ( call_number == 1 ) then
      write(tsh_file,*) " "
      write(tsh_file,"(1X,A,I10)") "Nuclear Step=", tsh_nucStep
      allocate(WFcis(ndets,all_states))
      call obtain_wavefunction(Xexc,WFcis,nstates,all_states,ndets,Ndim,NCO,Nvirt)

      ! Obtain Sigma Vectors now
      allocate(sigma(all_states,all_states))
      call obtain_sigma(WFcis,WFcis_old,C_scf,C_scf_old,sigma, &
                        ndets,all_states,M,NCO,Nvirt)
      WFcis_old = WFcis
      Sovl_old  = Sovl_now
      a_old = a; c_old = c; r_old = r
      C_scf_old = C_scf; sigma_old = sigma
   endif

   ! Obtain Forces
   call obtain_forces(Xexc,Eexc,C_scf,E_scf,M,Mlr,Nvirt,NCOlr,Ndim,nstates)

   ! Free memory of exchange basis 
   call basis_deinitLR()

   ! Step Nuclear
   tsh_nucStep = tsh_nucStep + 1

end subroutine TSHmain
