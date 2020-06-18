! Gets the LJS energy, which is the same as the conventional
! Lennard-Jones potential (eps is already 4*ε):
!
! E = eps * [ (sig/r)^12 - (sig/r)^6 ]
subroutine ljs_get_energy(energy)
   use LJ_switch_data, only: lj_atoms, mm_atoms, mmlj_sig, mmlj_eps

   implicit none
   LIODBLE, intent(out) :: energy

   integer :: iatom, jatom
   LIODBLE :: rterm, epsil

   energy = 0.0D0
   do iatom = 1, size(lj_atoms,1)
   do jatom = 1, size(mm_atoms,1)
      rterm = 0.5D0 * (mmlj_sig(mm_atoms(jatom)%mmtype) + &
                       lj_atoms(iatom)%sig)
      rterm = ( rterm / mm_atoms(jatom)%dist(iatom) ) ** 6

      epsil = sqrt (mmlj_eps(mm_atoms(jatom)%mmtype) * &
                    lj_atoms(iatom)%eps)
      
      ! eps is already stored as 4 * eps
      energy = energy + epsil * rterm * (rterm - 1.0D0)
   enddo
   enddo
end subroutine ljs_get_energy

! Calculates the term dE_LJ/dQ, with Q being the mulliken charge of
! the relevant QM atom (iatom). Also returns energy to avoid the above
! subroutine.
!
! The derivation is as follows, with deps being dε/dQ and dsig
! being dσ/dQ:
!
! dE_LJ/dQ = [ (sig/r)^12 - (sig/r)^6 ] deps + 
!            eps * [ 12 * sig^11/r^12 - 6 sig^5/r^6 ] * dsig
!
! And taking  x = (sig/r)^6 as a common factor:
!
! dE_dLJ/dQ = x * [ ( x - 1 ) * deps + (12 * x - 6) * eps * dsig / sig ]
! 
! Always take into acount that the ε stored in eps is 4*ε, and therefore
! deps = 4*deps.
! In addition, since sig = (sig_QM + sig_MM) / 2, dsig = dsig_QM/2; and
! eps = sqrt(eps_QM * eps_MM), deps = sqrt(eps_MM / eps_QM) * deps_QM / 2
subroutine ljs_get_dEdQ(energy, dE_dQ, iatom)
   use LJ_switch_data, only: lj_atoms, mm_atoms, mmlj_sig, mmlj_eps
   implicit none

   integer, intent(in)  :: iatom
   LIODBLE, intent(out) :: energy
   LIODBLE, intent(out) :: dE_dQ(:)
   
   integer :: jatom
   LIODBLE :: sigm, rterm, epsil, term_deps, term_dsig

   energy = 0.0D0
   dE_dQ  = 0.0D0
   do jatom = 1, size(mm_atoms,1)
      sigm  = 0.5D0 * (mmlj_sig(mm_atoms(jatom)%mmtype) + &
                       lj_atoms(iatom)%sig)
      rterm = ( sigm / mm_atoms(jatom)%dist(iatom) ) ** 6

      epsil = sqrt (mmlj_eps(mm_atoms(jatom)%mmtype) * &
                    lj_atoms(iatom)%eps)

      term_deps = (rterm - 1.0D0) * lj_atoms(iatom)%deps * epsil &
                  / lj_atoms(iatom)%eps
                 
      term_dsig = (12.0D0 * rterm - 6.0D0) * epsil &
                  * lj_atoms(iatom)%dsig / sigm

      ! This 0.5D0 extra comes from deps and dsig derivatives (see above)
      dE_dQ(iatom) = dE_dQ(iatom) + (term_deps + term_dsig) * rterm * 0.5D0

      ! Also calculates energy.
      energy = energy + epsil * rterm * (rterm - 1.0D0)
   enddo
   
end subroutine ljs_get_dEdQ