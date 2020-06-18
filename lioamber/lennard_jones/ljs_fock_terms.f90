subroutine ljs_add_fock_terms(fock, energ, rho, S_matrix, n_of_func, &
   atom_Z)
   use LJ_switch_data, only: n_lj_atoms, lj_atoms
   implicit none
   integer, intent(in)    :: n_of_func(:), atom_Z(:)
   LIODBLE, intent(in)    :: S_matrix(:,:), rho(:,:)
   LIODBLE, intent(inout) :: fock(:,:), energ

   integer :: ifunc, jfunc, iatom, f_idx
   LIODBLE :: atom_Q

   if (n_lj_atoms < 1) return

   do iatom = 1, n_lj_atoms

   ! Calculates the Mulliken charge on atom iatom
   do ifunc = 1, size(lj_atoms(iatom)%basis_id,1)
   f_idx = lj_atoms(iatom)%basis_id(ifunc)
   do jfunc = 1, size(fock,1) 
   atom_Q = atom_Q + rho(f_idx, jfunc) * S_matrix(f_idx, jfunc)
   enddo
   enddo
   atom_Q = lj_atoms(iatom)%Z - atom_Q 

   call lj_atoms(iatom)%set_eps_sig( atom_Q )


   enddo
end subroutine ljs_add_fock_terms
