!! VER EN AMBER: Opnq_LJ_atom_pair en SQM
! En params (nttyp, cn1, cn2)
! nttyp = ntypes*(ntypes+1)/2  (number of LJ types)
! ntypes es el total de combinaciones de LJ, por lo que
! solo se necesitan los elementos diagonales (i=i+1) de
! cn1 (A) y cn2 (B) de modo de sacar epsilon y sigma para
! cada átomo.

! En qm module, qmmm_struct: (iqmatoms, qm_mmpairs)
! iqmatoms tiene los indices de los átomos. 
! qm_mm_pairs tiene todos los pares mm-qm


!! Los MMTYPE van en las listas de epsilon y sigma. Estos vienen de
! ix, que es un array en memory_module.F90 aunque en ese mismo modulo
! tambien está atom type index en ese mismo modulo (importarlo)
! qmType=qmmm_struct%qm_atom_type(iqm)
! mmtype_for_iqm=qmmm_opnq%MM_atomType( qmmm_struct%iqmatoms(iqm) )
! jmm_index=qmmm_struct%qm_mm_pair_list(jmm)
! mmtype=qmmm_opnq%MM_atomType(jmm_index)


!! Number of pairs per QM atom. - length of pair_list.  
  ! integer :: qm_mm_pairs

  !! Non bond pair list for each QM atom
  ! integer, dimension(:), pointer :: qm_mm_pair_list => null()

#include "datatypes/datatypes.fh"
module LJ_switch_data
   implicit none

   type lj_atom
      integer :: idx    = 0
      integer :: Z      = 0
      integer :: mmtype = 0
      LIODBLE :: q1     = 0.0D0
      LIODBLE :: q2     = 0.0D0
      LIODBLE :: s1     = 0.0D0
      LIODBLE :: s2     = 0.0D0
      LIODBLE :: e1     = 0.0D0
      LIODBLE :: e2     = 0.0D0
      LIODBLE :: eps    = 0.0D0
      LIODBLE :: sig    = 0.0D0
      LIODBLE :: deps   = 0.0D0
      LIODBLE :: dsig   = 0.0D0
      integer, allocatable :: basis_id(:)

      contains
         procedure, pass :: set_eps_sig
         procedure, pass :: kill => destroy_lj
   end type lj_atom

   type mm_atom
      integer :: mmtype = 0
      LIODBLE, allocatable :: dist(:)

      contains
         procedure, pass :: kill => destroy_mm
   end type mm_atom

   LIODBLE :: k_fermi = 10.0D0
   integer :: n_lj_atoms
   type(lj_atom), allocatable :: lj_atoms(:)

   LIODBLE, allocatable       :: mmlj_eps(:) 
   LIODBLE, allocatable       :: mmlj_sig(:) 
   type(mm_atom), allocatable :: mm_atoms(:)

contains

   subroutine set_eps_sig(this, chrg)
      implicit none
      LIODBLE       , intent(in)    :: chrg
      class(lj_atom), intent(inout) :: this

      LIODBLE :: exp_term, inv_exp

      exp_term = exp( -k_fermi * ( chrg - 0.5D0 * (this%q2 + this%q1) ) )
      inv_exp  = 1.0D0 / (1.0D0 + exp_term)

      ! Base values
      this%sig = this%s1 + inv_exp * (this%s2 - this%s1)
      this%eps = this%e1 + inv_exp * (this%e2 - this%e1)

      ! Derivatives dSig/dQ and dEps/dQ
      inv_exp = k_fermi * exp_term * inv_exp * inv_exp

      this%dsig = inv_exp * (this%s2 - this%s1)
      this%deps = inv_exp * (this%e2 - this%e1)      

   end subroutine set_eps_sig

   subroutine destroy_lj(this)
      implicit none
      class(lj_atom), intent(inout) :: this
      
      if (allocated(this%basis_id)) deallocate(this%basis_id)
   end subroutine destroy_lj

   subroutine destroy_mm(this)
      implicit none
      class(mm_atom), intent(inout) :: this
      
      if (allocated(this%dist)) deallocate(this%dist)
   end subroutine destroy_mm

end module LJ_switch_data

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

subroutine ljs_substract_mm(energy, forces_qm, forces_mm)
   use LJ_switch_data, only: lj_atoms, mm_atoms, mmlj_sig, mmlj_eps
   use garcha_mod    , only: pos => r

   implicit none
   LIODBLE, intent(out) :: forces_qm(:)
   LIODBLE, intent(out) :: forces_mm(:)
   LIODBLE, intent(out) :: energy

   integer :: iatom, jatom
   LIODBLE :: rterm, rterm6, epsil

   energy = 0.0D0
   do iatom = 1, size(lj_atoms,1)
   do jatom = 1, size(mm_atoms,1)
      rterm = 0.5D0 * (mmlj_sig(mm_atoms(jatom)%mmtype) + &
                       mmlj_sig(lj_atoms(iatom)%mmtype))
      rterm6 = ( rterm / mm_atoms(jatom)%dist(iatom) ) ** 6

      epsil = sqrt (mmlj_eps(mm_atoms(jatom)%mmtype) * &
                    mmlj_eps(lj_atoms(iatom)%mmtype))
      
      ! eps is already stored as 4 * eps
      energy = energy + epsil * rterm6 * (rterm6 - 1.0D0)
   enddo
   enddo

end subroutine ljs_substract_mm


end module LJ_switch