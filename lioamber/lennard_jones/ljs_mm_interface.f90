subroutine ljs_set_params(eps_in, sig_in)
   use LJ_switch_data, only: mmlj_eps, mmlj_sig

   implicit none
   LIODBLE, intent(in) :: eps_in(:)
   LIODBLE, intent(in) :: sig_in(:)

   integer :: ntypes, itype

   if (allocated(mmlj_eps)) deallocate(mmlj_eps)
   if (allocated(mmlj_sig)) deallocate(mmlj_sig)

   ntypes = size(eps_in,1)
   allocate(mmlj_eps(ntypes), mmlj_sig(ntypes))
   
   do itype = 1, ntypes
      mmlj_eps(itype) = eps_in(itype)
      mmlj_sig(itype) = sig_in(itype)
   enddo
end subroutine

subroutine ljs_settle_mm(qm_types, mm_types, pos_qm, pos_mm)
   use LJ_switch_data, only: lj_atoms, mm_atoms
   implicit none
   integer, intent(in) :: qm_types(:)
   integer, intent(in) :: mm_types(:)
   LIODBLE, intent(in) :: pos_qm(:,:)
   LIODBLE, intent(in) :: pos_mm(:,:)

   integer :: iatom, jatom, n_solv, n_qm
   LIODBLE :: dist

   n_qm   = size(qm_types,1)
   n_solv = size(mm_types,1)
   if (allocated(mm_atoms)) then
      do iatom = 1, size(mm_atoms,1)
         call mm_atoms(iatom)%kill()
      enddo      
      deallocate(mm_atoms)
   endif
   allocate(mm_atoms(n_solv))

   do iatom = 1, size(lj_atoms,1)
      lj_atoms(iatom)%mmtype = qm_types(lj_atoms(iatom)%idx)
   enddo

   do iatom = 1, n_solv
      mm_atoms(iatom)%mmtype = mm_types(iatom)
      allocate(mm_atoms(iatom)%dist(size(lj_atoms,1)))

      do jatom = 1, size(lj_atoms,1)
         dist = ( pos_mm(1,iatom) - pos_qm(1,lj_atoms(jatom)%idx) ) * &
                ( pos_mm(1,iatom) - pos_qm(1,lj_atoms(jatom)%idx) )
         dist = ( pos_mm(2,iatom) - pos_qm(2,lj_atoms(jatom)%idx) ) * &
                ( pos_mm(2,iatom) - pos_qm(2,lj_atoms(jatom)%idx) ) + dist
         dist = ( pos_mm(3,iatom) - pos_qm(3,lj_atoms(jatom)%idx) ) * &
                ( pos_mm(3,iatom) - pos_qm(3,lj_atoms(jatom)%idx) ) + dist
         dist = sqrt(dist)

         mm_atoms(iatom)%dist(jatom) = dist
      enddo     
   enddo
end subroutine

! This subroutine substracts classical LJ terms. As such, all
! energy and gradient contributions are NEGATIVE.
subroutine ljs_substract_mm(energy, grads_qm, grads_mm, pos_qm, pos_mm)
   use LJ_switch_data, only: lj_atoms, mm_atoms, mmlj_sig, mmlj_eps

   implicit none
   LIODBLE, intent(in)  :: pos_qm(:,:)
   LIODBLE, intent(in)  :: pos_mm(:,:)
   LIODBLE, intent(out) :: grads_qm(:,:)
   LIODBLE, intent(out) :: grads_mm(:,:)
   LIODBLE, intent(out) :: energy

   integer :: iatom, jatom
   LIODBLE :: rterm, epsil, dE_dR_R, dx, dy, dz, dist

   do iatom = 1, size(lj_atoms,1)
   do jatom = 1, size(mm_atoms,1)
      ! Checks ε in order to avoid potential problems with
      ! σ = 0.0 (as in hydrogen atoms).
      if ((mmlj_eps(mm_atoms(jatom)%mmtype) > 0.0) .and. &
          (mmlj_eps(lj_atoms(iatom)%mmtype) > 0.0)) then
         dist  = mm_atoms(jatom)%dist(iatom)
         rterm = 0.5D0 * (mmlj_sig(mm_atoms(jatom)%mmtype) + &
                        mmlj_sig(lj_atoms(iatom)%mmtype))
         rterm = ( rterm / dist ) ** 6

         epsil = sqrt (mmlj_eps(mm_atoms(jatom)%mmtype) * &
                     mmlj_eps(lj_atoms(iatom)%mmtype))
         
         
         ! LJ epsilon is already stored as 4 * eps.
         ! Energy is SUBSTRACTED.
         energy = energy - epsil * rterm * (rterm - 1.0D0)

         ! Substracts classical LJ gradient terms. The first term is dE/dR
         ! divided by R, since it is a common factor between x, y and z.
         ! The negative sign is already included in the parenthesis, since
         ! the proper derivative is 6 - 12 * rterm.
         dE_dR_R = epsil * rterm * (12.0D0 * rterm - 6.0D0) / (dist * dist)

         dx = dE_dR_R * ( pos_qm(1,lj_atoms(iatom)%idx) - pos_mm(1,jatom) )
         dy = dE_dR_R * ( pos_qm(2,lj_atoms(iatom)%idx) - pos_mm(2,jatom) )
         dz = dE_dR_R * ( pos_qm(3,lj_atoms(iatom)%idx) - pos_mm(3,jatom) )

         grads_qm(1,lj_atoms(iatom)%idx) = grads_qm(1,lj_atoms(iatom)%idx) + dx
         grads_qm(2,lj_atoms(iatom)%idx) = grads_qm(2,lj_atoms(iatom)%idx) + dy
         grads_qm(3,lj_atoms(iatom)%idx) = grads_qm(3,lj_atoms(iatom)%idx) + dz

         grads_mm(1,jatom) = grads_mm(1,jatom) - dx
         grads_mm(2,jatom) = grads_mm(2,jatom) - dy
         grads_mm(3,jatom) = grads_mm(3,jatom) - dz
      endif
   enddo
   enddo
end subroutine ljs_substract_mm

subroutine ljs_gradients_qmmm(grads_qm, grads_mm, pos, n_qm)
   use LJ_switch_data, only: lj_atoms, mm_atoms, mmlj_sig, mmlj_eps

   implicit none
   integer, intent(in)  :: n_qm
   LIODBLE, intent(in)  :: pos(:,:)
   LIODBLE, intent(out) :: grads_qm(:,:)
   LIODBLE, intent(out) :: grads_mm(:,:)

   integer :: iatom, jatom
   LIODBLE :: rterm, epsil, dE_dR_R, dx, dy, dz, dist

   do iatom = 1, size(lj_atoms,1)
   do jatom = 1, size(mm_atoms,1)
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
   enddo
   enddo
end subroutine ljs_gradients_qmmm