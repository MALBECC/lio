! Adds LJS terms to Fock matrix, which are:
!
! dE_dRho_ij = Σa[ dE_dQa * dQa_dRho_ij ]
!            = Σa[- dE_dQa * S_ija ]
! 
! The S_ija terms are those terms of the overlap matrix
! which correspond to basis functions centered on atom a.
subroutine ljs_add_fock_terms(fock, energ, rho, S_matrix)
   use LJ_switch_data, only: n_lj_atoms, lj_atoms
   implicit none
   LIODBLE, intent(in)    :: S_matrix(:,:)
   LIODBLE, intent(in)    :: rho(:,:)
   LIODBLE, intent(inout) :: fock(:,:)
   LIODBLE, intent(inout) :: energ

   integer :: ifunc, jfunc, iatom, f_idx
   LIODBLE :: atom_Q, ener_atom, dEdQ

   if (n_lj_atoms < 1) return

   energ = 0.0D0
   do iatom = 1, n_lj_atoms

      ! Calculates the Mulliken charge on atom iatom
      atom_Q = 0.0D0
      do ifunc = 1, size(lj_atoms(iatom)%basis_id,1)
         f_idx = lj_atoms(iatom)%basis_id(ifunc)

         do jfunc = 1, size(S_matrix,1) 
            atom_Q = atom_Q + rho(f_idx, jfunc) * S_matrix(f_idx, jfunc)
         enddo
      enddo
      atom_Q = lj_atoms(iatom)%Z - atom_Q 

      ! Uses atomic charge to set the new LJ values for atom iatom, and
      ! calculates energy terms.
      call lj_atoms(iatom)%set_eps_sig( atom_Q )
      call ljs_get_dEdQ(ener_atom, dEdQ, iatom)
      energ = energ + ener_atom

      ! Calculates fock contributions.
      do ifunc = 1, size(lj_atoms(iatom)%basis_id,1)
         f_idx = lj_atoms(iatom)%basis_id(ifunc)

         do jfunc = 1, size(fock,1) 
            fock(f_idx, jfunc) = fock(f_idx, jfunc) - &
                                 S_matrix(f_idx, jfunc) * dEdQ
            if (f_idx /= jfunc) fock(jfunc, f_idx) = fock(jfunc, f_idx) - dEdQ &
                                                   * S_matrix(jfunc, f_idx)
         enddo
      enddo
   enddo
end subroutine ljs_add_fock_terms

subroutine ljs_add_fock_terms_op(fock, energ, rho, S_matrix, fockb, rhob)
   use LJ_switch_data, only: n_lj_atoms, lj_atoms
   implicit none
   LIODBLE, intent(in)    :: S_matrix(:,:)
   LIODBLE, intent(in)    :: rho(:,:)
   LIODBLE, intent(inout) :: fock(:,:)
   LIODBLE, intent(inout) :: energ

   LIODBLE, intent(in)   , optional :: rhob(:,:)
   LIODBLE, intent(inout), optional :: fockb(:,:)

   integer :: ifunc, jfunc, iatom, f_idx
   LIODBLE :: atom_Q, ener_atom, dEdQ

   if (n_lj_atoms < 1) return

   energ = 0.0D0
   do iatom = 1, n_lj_atoms

      ! Calculates the Mulliken charge on atom iatom
      atom_Q = 0.0D0
      do ifunc = 1, size(lj_atoms(iatom)%basis_id,1)
         f_idx = lj_atoms(iatom)%basis_id(ifunc)

         do jfunc = 1, size(S_matrix,1) 
            atom_Q = atom_Q + (rho(f_idx, jfunc) + rhob(f_idx, jfunc)) &
                            * S_matrix(f_idx, jfunc)
         enddo
      enddo
      atom_Q = lj_atoms(iatom)%Z - atom_Q 

      ! Uses atomic charge to set the new LJ values for atom iatom, and
      ! calculates energy terms.
      call lj_atoms(iatom)%set_eps_sig( atom_Q )
      call ljs_get_dEdQ(ener_atom, dEdQ, iatom)
      energ = energ + ener_atom

      ! Calculates fock contributions.
      do ifunc = 1, size(lj_atoms(iatom)%basis_id,1)
         f_idx = lj_atoms(iatom)%basis_id(ifunc)

         do jfunc = 1, size(fock,1) 
            fock(f_idx, jfunc) = fock(f_idx, jfunc) - &
                                 S_matrix(f_idx, jfunc) * dEdQ
            fockb(f_idx, jfunc) = fockb(f_idx, jfunc) - &
                                  S_matrix(f_idx, jfunc) * dEdQ
            if (f_idx /= jfunc) then
               fock(jfunc, f_idx) = fock(jfunc, f_idx) - dEdQ &
                                                       * S_matrix(jfunc, f_idx)
               fockb(jfunc, f_idx) = fockb(jfunc, f_idx) - dEdQ &
                                                       * S_matrix(jfunc, f_idx)
            endif
         enddo
      enddo
   enddo
end subroutine ljs_add_fock_terms_op