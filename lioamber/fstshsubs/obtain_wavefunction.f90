subroutine obtain_wavefunction(lrcoef,wfunc,nstates,all_states,ndets,Ndim,NCO,Nvirt)
! lrcoef = Excited States Coefficients obtained from LR
! wfunc  = Wavefunction in CIS format
! nstates= only excited states
! all_states = excited + ground states
! ndets  = number of determinants, including ground states
! Ndim   = Dimension of LR coeff.
! NCO    = number of molecular occupied orbitals
! Nvirt  = number of molecular unoccupied orbitals

! ##################### EQUATION #########################
! \Phi^I = det_gs + \sum_ia{C^I_ia * det_ia}

! \Phi^I = CIS wavefunction of I state
! det_gs = ground state determinant
! i      = molecular occupied orbital
! a      = molecular unoccupied orbital
! C^I_ia = LR coef of I state, element i-a
! det_ia = Excited state determinant, switch i by a
! ########################################################
   implicit none

   integer, intent(in)  :: nstates, all_states, ndets, Ndim, NCO, Nvirt
   LIODBLE, intent(in)  :: lrcoef(Ndim,nstates)
   LIODBLE, intent(out) :: wfunc(ndets,all_states)

   LIODBLE :: norm, val
   integer :: ii, iocc, jvirt, ind, ind2, NCOc

   norm = 1.0d0/dsqrt(2.0d0)
   NCOc = NCO - 1

   ! Ground States
   Wfunc(1,1) = 1.0d0
   Wfunc(2:ndets,1) = 0.0d0

   ! Excited States
   do ii=2,all_states
      wfunc(1,ii) = 0.0d0 ! Contribution of GS to ES

      do iocc=0,NCOc
      do jvirt=1,Nvirt
         ind  = iocc*Nvirt+jvirt
         ind2 = (NCOc-iocc)*Nvirt+jvirt-1
         val  = lrcoef(ind,ii-1) * norm
         wfunc(1+2*ind2+1,ii) = val
         wfunc(1+2*ind2+2,ii) = val
      enddo
      enddo
   enddo
 
! TODO: maybe we can orthonormalize the CIS wavefunction

end subroutine obtain_wavefunction
