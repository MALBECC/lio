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
                       ljatoms(iatom)%sig)
      rterm = ( rterm / mm_atoms(jatom)%dist(iatom) ) ** 6

      epsil = sqrt (mmlj_eps(mm_atoms(jatom)%mmtype) * &
                    ljatoms(iatom)%eps)
      
      ! eps is already stored as 4 * eps
      energy = energy + epsil * rterm * (rterm - 1.0D0)
   enddo
   enddo
end subroutine ljs_get_energy


subroutine ljs_get_dEdQ()
   implicit none

end subroutine ljs_get_dEdQ