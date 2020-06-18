subroutine ljs_settle_mm(qm_types, mm_types)
   use LJ_switch_data, only: n_lj_atoms, lj_atoms, mm_atoms
   use garcha_mod,     only: r, natoms
   implicit none
   integer, intent(in) :: qm_types(:)
   integer, intent(in) :: mm_types(:)

   integer :: iatom, jatom, n_solv
   LIODBLE :: dist

   n_solv = size(mm_types,1)
   if (allocated(mm_atoms)) then
      do iatom = 1, size(mm_atoms,1)
         call mm_atoms(iatom)%kill()
      enddo      
      deallocate(mm_atoms)
   endif
   allocate(mm_atoms(n_solv))

   do iatom = 1, n_lj_atoms
      lj_atoms(iatom)%mmtype = qm_types(lj_atoms(iatom)%idx)
   enddo

   do iatom = 1, n_solv
      mm_atoms(iatom)%mmtype = mm_types(iatom)
      allocate(mm_atoms(iatom)%dist(n_lj_atoms))

      do jatom = 1, n_lj_atoms
         dist = (r(natom + iatom,1) - r(lj_atoms(iatom)%idx,1)) * &
                (r(natom + iatom,1) - r(lj_atoms(iatom)%idx,1))
         dist = (r(natom + iatom,2) - r(lj_atoms(iatom)%idx,2)) * &
                (r(natom + iatom,2) - r(lj_atoms(iatom)%idx,2)) + dist
         dist = (r(natom + iatom,3) - r(lj_atoms(iatom)%idx,3)) * &
                (r(natom + iatom,3) - r(lj_atoms(iatom)%idx,3)) + dist
         dist = sqrt(dist)

         mm_atoms(iatom)%dist(jatom) = dist
      enddo     
   enddo
end subroutine

subroutine ljs_substract_mm(energy, grads_qm, grads_mm)
   use LJ_switch_data, only: lj_atoms, mm_atoms, mmlj_sig, mmlj_eps
   use garcha_mod    , only: pos => r, natom

   implicit none
   LIODBLE, intent(out) :: grads_qm(:,3)
   LIODBLE, intent(out) :: grads_mm(:,3)
   LIODBLE, intent(out) :: energy

   integer :: iatom, jatom
   LIODBLE :: rterm, epsil, dE_dR_R, dx, dy, dz

   energy = 0.0D0
   do iatom = 1, size(lj_atoms,1)
   do jatom = 1, size(mm_atoms,1)
      dist  = mm_atoms(jatom)%dist(iatom)
      rterm = 0.5D0 * (mmlj_sig(mm_atoms(jatom)%mmtype) + &
                       mmlj_sig(lj_atoms(iatom)%mmtype))
      rterm = ( rterm / dist ) ** 6

      epsil = sqrt (mmlj_eps(mm_atoms(jatom)%mmtype) * &
                    mmlj_eps(lj_atoms(iatom)%mmtype))
      
      ! LJ epsilon is already stored as 4 * eps
      energy = energy + epsil * rterm * (rterm - 1.0D0)

      ! Substracts classical LJ gradient terms. The first term is dE/dR
      ! divided by R, since it is a common factor between x, y and z.
      dE_dR_R = epsil * rterm * (6.0D0 - 12.0D0 * rterm) / (dist * dist)


      dx = dE_dR_R * (r(lj_atoms(iatom)%idx,1) - r(natom + jatom,1))
      dy = dE_dR_R * (r(lj_atoms(iatom)%idx,2) - r(natom + jatom,2))
      dz = dE_dR_R * (r(lj_atoms(iatom)%idx,3) - r(natom + jatom,3))

      grads_qm(1,lj_atoms(iatom)%idx) = grads_qm(1,lj_atoms(iatom)%idx) + dx
      grads_qm(2,lj_atoms(iatom)%idx) = grads_qm(2,lj_atoms(iatom)%idx) + dy
      grads_qm(3,lj_atoms(iatom)%idx) = grads_qm(3,lj_atoms(iatom)%idx) + dz

      grads_mm(1,jatom) = grads_mm(1,jatom) - dx
      grads_mm(2,jatom) = grads_mm(2,jatom) - dy
      grads_mm(3,jatom) = grads_mm(3,jatom) - dz
   enddo
   enddo
end subroutine ljs_substract_mm
