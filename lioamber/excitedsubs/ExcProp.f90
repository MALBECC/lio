subroutine ExcProp(CoefA,CoefB,EneA,EneB,Etot)
! This is the main routine to calculate excited states properties
! For now this perform:
! - Linear Response Calculation


! Inputs
! - CoefA: Molecular Orbitals Coefficient of alpha
! - CoefB: Molecular Orbitals COefficient of beta
! - EneA: Molecular Orbitals Energy of alpha
! - EneB: Molecular Orbitals Energy of beta
use garcha_mod, only: OPEN, NCO, PBE0
use excited_data, only: lresp, nstates, libint_recalc, fittExcited
use basis_data, only: M, c_raw
   implicit none

   double precision, intent(in) :: CoefA(:,:), CoefB(:,:)
   double precision, intent(in) :: EneA(:), EneB(:)
   double precision, intent(inout) :: Etot

   integer :: NCOlr, Mlr, Nvirt, Ndim
   double precision, allocatable :: C_scf(:,:), E_scf(:)
   double precision, allocatable :: Xexc(:,:), Eexc(:)

   if (lresp .eqv. .false.) return
   if (OPEN  .eqv. .true. ) then 
      print*, "Linear Response doesn't work in Open shell"
      stop
   endif

   if ( (.not. PBE0) .and.  (.not. fittExcited) ) then
      call g2g_timer_sum_start('Libint init')
      call g2g_libint_init(c_raw,libint_recalc)
      call g2g_timer_sum_stop('Libint init')
   endif

   ! This routine applies the FCA method
   ! NCO = number of occupied molecular orbitals
   ! Nvirt = number of virtual molecular orbitals
   ! Ndim = dimension of Excited Matrix = (NCOxNvirt)^2
   ! C_scf, E_scf = Molecular Orbital Coeff. and Energy
   call fcaApp(CoefA,EneA,C_scf,E_scf,NCO,M,NCOlr,Mlr,Nvirt,Ndim)

   ! This routine form matrices for change basis
   call basis_initLR(C_scf,M,Mlr,NCOlr,Nvirt)

   ! Save density and derivatives values of Ground State
   call g2g_timer_start("Save GS Density")
   call g2g_saverho( )
   call g2g_timer_stop("Save GS Density")

   ! Linear Response Calculation
   ! This routine obtain the Excitation Energy and
   ! Transition Vectors
   allocate(Xexc(Ndim,nstates),Eexc(nstates))
   call g2g_timer_start("Linear Response")
   call linear_response(C_scf,E_scf,Xexc,Eexc,M,Mlr,Nvirt,NCOlr,Ndim,0)
   call g2g_timer_stop("Linear Response")

   ! Relaxed Density Matrix of one Excited State
   call RelaxedDensity(Xexc,Eexc,C_scf,E_scf,M,Nvirt,NCO,Ndim,nstates)






   ! Deinitialization
   call basis_deinitLR()
  
end subroutine ExcProp
