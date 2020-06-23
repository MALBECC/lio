
subroutine ljs_gradients_qmmm(grads_qm, grads_mm, pos, n_qm)
   use LJ_switch_data, only: lj_atoms, mm_atoms, mmlj_sig, mmlj_eps, n_lj_atoms

   implicit none
   integer, intent(in)  :: n_qm
   LIODBLE, intent(in)  :: pos(:,:)
   LIODBLE, intent(out) :: grads_qm(:,:)
   LIODBLE, intent(out) :: grads_mm(:,:)

   integer :: iatom, jatom
   LIODBLE :: rterm, epsil, dE_dR_R, dx, dy, dz, dist

   if (n_lj_atoms < 1) return
   do iatom = 1, size(lj_atoms,1)
   do jatom = 1, size(mm_atoms,1)
      if ((abs(mmlj_eps(mm_atoms(jatom)%mmtype)) > 0.0) .and. &
          (abs(lj_atoms(iatom)%eps) > 0.0)) then
         dist  = mm_atoms(jatom)%dist(iatom)
         rterm = 0.5D0 * (mmlj_sig(mm_atoms(jatom)%mmtype) + &
                        lj_atoms(iatom)%sig)
         rterm = ( rterm / dist ) ** 6

         epsil = sqrt (mmlj_eps(mm_atoms(jatom)%mmtype) * &
                     lj_atoms(iatom)%eps)
         
         ! Substracts classical LJ gradient terms. The first term is dE/dR
         ! divided by R, since it is a common factor between x, y and z.
         dE_dR_R = epsil * rterm * (6.0D0 - 12.0D0 * rterm) / (dist * dist)


         dx = dE_dR_R * (pos(lj_atoms(iatom)%idx,1) - pos(n_qm + jatom,1))
         dy = dE_dR_R * (pos(lj_atoms(iatom)%idx,2) - pos(n_qm + jatom,2))
         dz = dE_dR_R * (pos(lj_atoms(iatom)%idx,3) - pos(n_qm + jatom,3))

         grads_qm(1,lj_atoms(iatom)%idx) = grads_qm(1,lj_atoms(iatom)%idx) + dx
         grads_qm(2,lj_atoms(iatom)%idx) = grads_qm(2,lj_atoms(iatom)%idx) + dy
         grads_qm(3,lj_atoms(iatom)%idx) = grads_qm(3,lj_atoms(iatom)%idx) + dz

         grads_mm(1,jatom) = grads_mm(1,jatom) - dx
         grads_mm(2,jatom) = grads_mm(2,jatom) - dy
         grads_mm(3,jatom) = grads_mm(3,jatom) - dz
      endif
   enddo
   enddo
end subroutine ljs_gradients_qmmm