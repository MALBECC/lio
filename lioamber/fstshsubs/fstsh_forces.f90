subroutine obtain_forces(Xexc,Eexc,C_scf,E_scf,M,Mlr,Nvirt,NCOlr,Ndim,nstates)
use excitedsubs , only: RelaxedDensity, forcesexc
use excited_data, only: root, excited_forces
use fstsh_data  , only: current_state
   implicit none

   integer, intent(in) :: M, Mlr, Nvirt, NCOlr, Ndim, nstates
   LIODBLE, intent(in) :: C_scf(M,Mlr), E_scf(M)
   LIODBLE, intent(in) :: Xexc(Ndim,nstates), Eexc(nstates)

   LIODBLE, dimension(:), allocatable   :: Zvec, Qvec
   LIODBLE, dimension(:,:), allocatable :: Gxc, rhoEXC, Pdif, Trans

   ! This routine only calculates the excited state force
   if ( current_state == 1 ) then
      root = 0
      excited_forces = .false.
      return
   endif

   ! Set variables needed to forces calculate
   root = current_state - 1
   excited_forces = .true.

   allocate(Zvec(Ndim),Qvec(Ndim),Gxc(M,M))
   allocate(rhoEXC(M,M),Pdif(M,M),Trans(M,M))
   call RelaxedDensity(Xexc,C_scf,E_scf,Zvec,Qvec,Gxc, &
                       rhoEXC,Pdif,Trans,M,Mlr,Nvirt,NCOlr,Ndim,nstates)

   ! Excited States Forces: This save forces in excited_data module
   call forcesexc(rhoEXC,Pdif,Zvec,Trans,Qvec,Gxc,Xexc,Eexc, &
                  C_scf,E_scf,M,Mlr,Ndim,NCOlr,nstates)

   deallocate(Zvec,Qvec,Gxc,rhoEXC,Pdif,Trans)
end subroutine obtain_forces
