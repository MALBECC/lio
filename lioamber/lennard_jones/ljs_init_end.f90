subroutine ljs_initialise(atom_Z, atom_of_func, is_doing_ljs)
   use LJ_switch_data, only: n_lj_atoms, lj_atoms
   
   implicit none
   integer, intent(in) :: atom_Z(:)
   integer, intent(in) :: atom_of_func(:)
   logical, intent(out), optional :: is_doing_ljs

   integer :: iatom, ifunc, f_count

   if (present(is_doing_ljs)) is_doing_ljs = .false.
   if (n_lj_atoms < 1) return

   if (present(is_doing_ljs)) is_doing_ljs = .true.
   do iatom = 1, n_lj_atoms
      lj_atoms(iatom)%Z = atom_Z(lj_atoms(iatom)%idx)

      f_count = 0
      do ifunc = 1, size(atom_of_func,1)
         if (atom_of_func(ifunc) == lj_atoms(iatom)%idx) then
            f_count = f_count +1
         endif
      enddo

      allocate(lj_atoms(iatom)%basis_id(f_count))
      f_count = 0
      do ifunc = 1, size(atom_of_func,1)
         if (atom_of_func(ifunc) == lj_atoms(iatom)%idx) then
            f_count = f_count +1
            lj_atoms(iatom)%basis_id(f_count) = ifunc
         endif
      enddo
   enddo
end subroutine ljs_initialise

subroutine ljs_finalise()
   use LJ_switch_data, only: mm_atoms, lj_atoms, mmlj_eps, mmlj_sig, n_lj_atoms
   implicit none
   integer :: iatom

   if (allocated(lj_atoms)) then
      do iatom = 1, size(lj_atoms,1)
         call lj_atoms(iatom)%kill()
      enddo 
      deallocate(lj_atoms)
   endif

   if (allocated(mm_atoms)) then
      do iatom = 1, size(mm_atoms,1)
         call mm_atoms(iatom)%kill()
      enddo
      deallocate(mm_atoms)
   endif

   if (allocated(mmlj_eps)) deallocate(mmlj_eps)
   if (allocated(mmlj_sig)) deallocate(mmlj_sig)
   n_lj_atoms = 0
end subroutine ljs_finalise
