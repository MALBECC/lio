subroutine ExcProp(Etot, CoefA, EneA, CoefB, EneB)
! This is the main routine to calculate excited states properties
! For now this perform:
! - Linear Response Calculation


! Inputs
! - CoefA: Molecular Orbitals Coefficient of alpha
! - CoefB: Molecular Orbitals COefficient of beta
! - EneA: Molecular Orbitals Energy of alpha
! - EneB: Molecular Orbitals Energy of beta
use garcha_mod  , only: OPEN, NCO, Pmat_vec
use excited_data, only: lresp, nstates, root, pack_dens_exc, second_LR, & 
                        Tdip_save, save_tlr, state_LR
use basis_data  , only: M
use td_data     , only: timedep
   implicit none

   LIODBLE, intent(in)           :: CoefA(:,:), EneA(:)
   LIODBLE, intent(in), optional :: CoefB(:,:), EneB(:)
   LIODBLE, intent(inout)        :: Etot

   integer :: NCOlr, Mlr, Nvirt, Ndim, ii
   LIODBLE, allocatable :: C_scf(:,:), E_scf(:)
   LIODBLE, allocatable :: Xexc(:,:), Eexc(:)
   LIODBLE, allocatable :: Zvec(:), Qvec(:), Gxc(:,:)
   LIODBLE, allocatable :: rhoEXC(:,:), Pdif(:,:), Trans(:,:)
   LIODBLE, allocatable :: Coef1(:,:), Xflr(:,:), Eflr(:)
   LIODBLE, allocatable :: Tdip0(:,:)

   if (lresp .eqv. .false.) return
   if (OPEN  .eqv. .true. ) then 
      print*, "Linear Response doesn't work in Open shell"
      stop
   endif

   ! DUMMY LINE FOR WARNINGS
   if (present(EneB) .and. present(CoefB)) print*, "Open shell not supported."

   ! Truncated MOs
   call truncated_MOs(CoefA,EneA,C_scf,E_scf,NCO,M,NCOlr,Mlr,Nvirt,Ndim)

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

   ! Saving Transition vectors:TODO: move this to other place
   if ( save_tlr ) then
      write(*,"(1X,A,I3)") "Saving state:", state_LR
      open(unit=456,file="vectors_root")
      write(456,*) Eexc(state_LR)
      do ii=1,Ndim
         write(456,*) Xexc(ii,state_LR)
      enddo
      close(456)
   endif

   ! Second Linear Response Calculation
   ! This routine obtain the Excited State Absorption
   if ( second_LR ) then
      ! Deinitialization of change basis
      call basis_deinitLR()
      allocate(Coef1(M,M),Xflr(Ndim,nstates),Eflr(nstates),Tdip0(nstates,3))

      ! For perturbed variables: MOS, density
      !   This put the perturbed density into Pmat_vec
      call get_perturbed(Coef1,CoefA,Xexc,M,Nvirt,NCO,Ndim,nstates)
 
      ! Truncated MOs, init change basis and density
      call truncated_MOs(Coef1,EneA,C_scf,E_scf,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
      call basis_initLR(C_scf,M,Mlr,NCOlr,Nvirt)
      call g2g_saverho( )

      ! Save first LR variables
      Xflr = Xexc; Eflr = Eexc; Tdip0 = Tdip_save
      deallocate(Xexc,Eexc); allocate(Xexc(Ndim,nstates),Eexc(nstates))
      Xexc = 0.0d0; Eexc = 0.0d0; Tdip_save = 0.0d0
    
      ! Second LR
      call linear_response(C_scf,E_scf,Xexc,Eexc,M,Mlr,Nvirt,NCOlr,Ndim,1) ! realizo el slr
      call fca_restored(CoefA,EneA,C_scf,E_scf,Xexc,M,Mlr,Nvirt,NCO,NCOlr,&
                        Ndim,nstates)
      
      ! Final excited spectra
      call excited_absorption(Xflr,Eflr,Xexc,Eexc,Coef1,Tdip0,M,NCO,Nvirt,Ndim,nstates) ! genero el spectro
      deallocate(Coef1,Xflr,Eflr,Tdip0)
   endif

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

   ! Cubegen of excited states properties
   call getcubegen_excited(Trans,Pdif,rhoExc,C_scf,M)

   ! Excited States Forces: This save forces in excited_data module
   call forcesexc(rhoEXC,Pdif,Zvec,Trans,Qvec,Gxc,Xexc,Eexc, &
                  C_scf,E_scf,M,Mlr,Ndim,NCOlr,nstates)

   ! Time Dependent Real Time with Excited State
   if ( timedep == 1 ) then
      print*, "RT-TDDFT with Excited Density State", root
      if (allocated(pack_dens_exc)) then
         Pmat_vec = pack_dens_exc
      else
         print*, "Excited State Density was not saved"
         stop
      endif
   endif

   ! Deinitialization and Free Memory
   call basis_deinitLR()
   deallocate(Zvec,Qvec,Gxc,rhoEXC,Pdif,Trans)
   deallocate(Xexc,Eexc)
  
end subroutine ExcProp
