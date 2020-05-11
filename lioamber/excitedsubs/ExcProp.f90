subroutine ExcProp(CoefA, CoefB, EneA, EneB, Etot)
! This is the main routine to calculate excited states properties
! For now this perform:
! - Linear Response Calculation


! Inputs
! - CoefA: Molecular Orbitals Coefficient of alpha
! - CoefB: Molecular Orbitals COefficient of beta
! - EneA: Molecular Orbitals Energy of alpha
! - EneB: Molecular Orbitals Energy of beta
use garcha_mod, only: OPEN, NCO
use excited_data, only: lresp, nstates, libint_recalc, fittExcited
use extern_functional_data, only: libint_inited
use extern_functional_subs, only: libint_init
use basis_data, only: M, c_raw
   implicit none

   LIODBLE, intent(in) :: CoefA(:,:), CoefB(:,:)
   LIODBLE, intent(in) :: EneA(:), EneB(:)
   LIODBLE, intent(inout) :: Etot

   integer :: NCOlr, Mlr, Nvirt, Ndim
   LIODBLE, allocatable :: C_scf(:,:), E_scf(:)
   LIODBLE, allocatable :: Xexc(:,:), Eexc(:)
   LIODBLE, allocatable :: Zvec(:), Qvec(:), Gxc(:,:)
   LIODBLE, allocatable :: rhoEXC(:,:), Pdif(:,:), Trans(:,:)

   if (lresp .eqv. .false.) return
   if (OPEN  .eqv. .true. ) then 
      print*, "Linear Response doesn't work in Open shell"
      stop
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

   ! This routine obtain the new indexes in order to delete FCA
   call fca_restored(CoefA,EneA,C_scf,E_scf,Xexc,M,Mlr,Nvirt,NCO,NCOlr,&
                     Ndim,nstates)

   ! This routine obtain Non-Adiabatic Coupling Vectors and
   ! evolution coefficients
   call tsh_probabilities(C_scf,E_scf,Xexc,Eexc,NCOlr,M,Mlr,Ndim,Nvirt,&
                          Etot,nstates)

   ! Relaxed Density Matrix of one Excited State
   allocate(Zvec(Ndim),Qvec(Ndim),Gxc(M,M))
   allocate(rhoEXC(M,M),Pdif(M,M),Trans(M,M))
   ! rhoEXC = Relaxed Density Matrix of Excited States root
   ! Pdif   = Difference Density Matrix
   call RelaxedDensity(Xexc,C_scf,E_scf,Zvec,Qvec,Gxc, &
                       rhoEXC,Pdif,Trans,M,Mlr,Nvirt,NCOlr,Ndim,nstates)

   ! Excited States Forces: This save forces in excited_data module
   call forcesexc(rhoEXC,Pdif,Zvec,Trans,Qvec,Gxc,Xexc,Eexc, &
                  C_scf,E_scf,M,Mlr,Ndim,NCOlr,nstates)

   





   ! Deinitialization and Free Memory
   call basis_deinitLR()
   deallocate(Zvec,Qvec,Gxc,rhoEXC,Pdif,Trans)
   deallocate(Xexc,Eexc)
  
end subroutine ExcProp
