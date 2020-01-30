subroutine ExcProp(CoefA,CoefB,EneA,EneB,Etot)
! This is the main routine to calculate excited states properties
! For now this perform:
! - Linear Response Calculation


! Inputs
! - CoefA: Molecular Orbitals Coefficient of alpha
! - CoefB: Molecular Orbitals COefficient of beta
! - EneA: Molecular Orbitals Energy of alpha
! - EneB: Molecular Orbitals Energy of beta
use garcha_mod, only: OPEN, NCO, PBE0, Pmat_vec
use excited_data, only: lresp, nstates, libint_recalc, fittExcited, &
                        excited_forces, pack_dens_exc
use basis_data, only: M, c_raw
   implicit none

   double precision, intent(in) :: CoefA(:,:), CoefB(:,:)
   double precision, intent(in) :: EneA(:), EneB(:)
   double precision, intent(inout) :: Etot

   integer :: NCOlr, Mlr, Nvirt, Ndim
   double precision, allocatable :: C_scf(:,:), E_scf(:)
   double precision, allocatable :: Xexc(:,:), Eexc(:)
   double precision, allocatable :: Zvec(:), Qvec(:), Gxc(:,:)
   double precision, allocatable :: rhoEXC(:,:), Pdif(:,:), Trans(:,:)

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
   allocate(Zvec(Ndim),Qvec(Ndim),Gxc(M,M))
   allocate(rhoEXC(M,M),Pdif(M,M),Trans(M,M))
   ! rhoEXC = Relaxed Density Matrix of Excited States root
   ! Pdif   = Difference Density Matrix
   call RelaxedDensity(Xexc,Eexc,C_scf,E_scf,Zvec,Qvec,Gxc, &
                       rhoEXC,Pdif,Trans,M,Mlr,Nvirt,NCOlr,Ndim,nstates)

   ! Excited States Forces: This save forces in excited_data module
   call forcesexc(rhoEXC,Pdif,Zvec,Trans,Qvec,Gxc,Xexc,Eexc, &
                  C_scf,E_scf,M,Mlr,Ndim,NCOlr,nstates)


   ! Check if we perform analysis with excited density matrix or not
   if ( excited_forces ) Pmat_vec = pack_dens_exc
   
   





   ! Deinitialization
   call basis_deinitLR()
  
end subroutine ExcProp
