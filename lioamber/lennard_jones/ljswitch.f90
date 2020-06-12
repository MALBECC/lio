#include "../datatypes/datatypes.fh"
module LJ_switch
   implicit none
   private
   public :: doing_ljs
   public :: ljs_input_read
   public :: ljs_add_fock_terms

contains

function doing_ljs() result(is_doing)
   use LJ_switch_data, only: n_lj_atoms
   
   implicit none
   logical :: is_doing

   is_doing = .false.
   if (n_lj_atoms > 0) is_doing = .true.
   return
end function doing_ljs

subroutine ljs_input_read(input_UID, verbose_lvl)
   use LJ_switch_data, only: n_lj_atoms, lj_atoms
   implicit none
   integer, intent(in)    :: input_UID
   integer, intent(in)    :: verbose_lvl
    
   character(len=10) :: buffer
   character(len=50) :: print_fmt
   integer           :: ios, iatom
    
   rewind(input_UID)
   ios = 0
   do while ((trim(buffer) /= "{LJSWITCH}") .and. (ios == 0) )
      read(input_UID,'(A10)', iostat = ios) buffer
   enddo
 
   ! If ios < 0, found EOF. No LJ Switch input provided.
   if (ios < 0) return
   write(*,'(A)') ""
   write(*,'(A)') "== LJ Switch =="

   ! Checks the number of input lines and rewinds for further read.
   iatom = -1
   do while ((trim(buffer) /= "{END}") .and. (ios == 0) )
      iatom = iatom + 1
      read(input_UID,'(A10)', iostat = ios) buffer
   enddo
   n_lj_atoms = iatom
   
   write(*,'(A30,I3)') "Input found. Number of atoms: ", n_lj_atoms
   rewind(input_UID)
   do while ((trim(buffer) /= "{LJSWITCH}") .and. (ios == 0) )
      read(input_UID,'(A10)', iostat = ios) buffer
   enddo

   if (allocated(lj_atoms)) deallocate(lj_atoms)
   allocate(lj_atoms(n_lj_atoms))

   ! Starts reading LJ switch data.
   do iatom = 1, n_lj_atoms
      read(input_UID,*) lj_atoms(iatom)%idx, lj_atoms(iatom)%q1, &
                        lj_atoms(iatom)%q2 , lj_atoms(iatom)%s1, &
                        lj_atoms(iatom)%s2 , lj_atoms(iatom)%e1, &
                        lj_atoms(iatom)%e2
      if (verbose_lvl > 2) then
         print_fmt = "(A5, 1x, I3, 5x, A6, I3)"
         write(*,print_fmt) "Atom: ", iatom, "Index: " lj_atoms(iatom)%idx

         print_fmt = "(A5, 1x, I1, A6, F12.6, A10, F12.6, A12, F12.6)"
         write(*,print_fmt) "Type", 1, " - Q: ", lj_atoms(iatom)%q1,         &
                            " - Sigma: ", lj_atoms(iatom)%s1," - Epsilon: ", &
                            lj_atoms(iatom)%e1
         write(*,print_fmt) "Type", 2, " - Q: ", lj_atoms(iatom)%q2,         &
                            " - Sigma: ", lj_atoms(iatom)%s2," - Epsilon: ", &
                            lj_atoms(iatom)%e2
      endif
   enddo
   write(*,'(A)') ""
end subroutine ljs_input_read

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


end module LJ_switch
