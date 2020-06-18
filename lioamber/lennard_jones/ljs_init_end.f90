subroutine ljs_initialise(eps_in, sig_in)
   use LJ_switch_data, only: n_lj_atoms, lj_atoms, mmlj_eps, mmlj_sig
   use garcha_mod    , only: atom_Z => Iz, atom_of_func => Nuc

   implicit none
   LIODBLE, intent(in) :: eps_in(:)
   LIODBLE, intent(in) :: sig_in(:)

   integer :: iatom, ifunc, f_count, ntypes, itype

   do iatom = 1, n_lj_atoms
      lj_atoms(iatom)%Z = atom_Z(lj_atoms(iatom)%idx)

      f_count = 0
      do ifunc = 1, size(atom_of_func,1)
         if (atom_of_func(iatom) == lj_atoms(iatom)%idx) then
            f_count = f_count +1
         endif
      enddo

      allocate(lj_atoms(iatom)%basis_id(f_count))
      f_count = 0
      do ifunc = 1, size(atom_of_func,1)
         if (atom_of_func(iatom) == lj_atoms(iatom)%idx) then
            f_count = f_count +1
            lj_atoms(iatom)%basis_id(f_count) = ifunc
         endif
      enddo
   enddo

   if (allocated(mmlj_eps)) deallocate(mmlj_eps)
   if (allocated(mmlj_sig)) deallocate(mmlj_sig)

   ntypes = size(eps_in,1)
   allocate(mmlj_eps(ntypes), mmlj_sig(ntypes))
   
   do itype = 1, ntypes
      mmlj_eps(itype) = eps_in(itype)
      mmlj_sig(itype) = sig_in(itype)
   enddo
end subroutine ljs_initialise