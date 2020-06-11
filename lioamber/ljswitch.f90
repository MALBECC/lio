!! VER EN AMBER
! En params:
 ! nttyp = ntypes*(ntypes+1)/2  (number of LJ type pairs)
 ! cn1 -> A, cn2 -> B

! En qm module:
  ! qmmm_struct%iqmatoms -> indice de los atomos QM

  !! Number of pairs per QM atom. - length of pair_list.  
  ! integer :: qm_mm_pairs

  !! Non bond pair list for each QM atom
  ! integer, dimension(:), pointer :: qm_mm_pair_list => null()

  !! atomic numbers of MM atoms included in QM-MM pairs (only used for PM3/MM*)
  !! qm_mm_pairs long, allocated in read_qmmm_nm_and_alloc for SQM external charges
  !! allocated in ??? for sander QM/MM
  ! integer, dimension(:), pointer :: qm_mm_pair_atom_numbers  => null()

#include "datatypes/datatypes.fh"
module LJ_switch_data
   implicit none

   type lj_atom
      integer :: idx  = 0
      integer :: Z    = 0
      LIODBLE :: q1   = 0.0D0
      LIODBLE :: q2   = 0.0D0
      LIODBLE :: s1   = 0.0D0
      LIODBLE :: s2   = 0.0D0
      LIODBLE :: e1   = 0.0D0
      LIODBLE :: e2   = 0.0D0
      LIODBLE :: eps  = 0.0D0
      LIODBLE :: sig  = 0.0D0
      LIODBLE :: deps = 0.0D0
      LIODBLE :: dsig = 0.0D0

      integer, allocatable :: basis_id(:)

      contains
         procedure, pass :: set_eps_sig
   end type lj_atom

   type mm_atom
      LIODBLE :: eps  = 0.0D0
      LIODBLE :: sig  = 0.0D0

      LIODBLE, allocatable :: dist(:)
   end type mm_atom

   integer :: n_lj_atoms
   logical :: ljs_initialised = .false.
   LIODBLE :: k_fermi = 10.0D0

   type(lj_atom), allocatable :: lj_atoms(:)
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

subroutine ljs_input_read(input_UID, verbose_lvl, do_mullik)
   use LJ_switch_data, only: n_lj_atoms, lj_atoms
   implicit none
   integer, intent(in)    :: input_UID
   integer, intent(in)    :: verbose_lvl
   logical, intent(inout) :: do_mullik
    
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
      if (verbose_lvl > 3) then
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

   if (n_lj_atoms > 0) do_mullik = .true.
   
end subroutine ljs_input_read


subroutine ljs_initialise(atom_of_func, atom_Z)
   use LJ_switch_data, only: n_lj_atoms, lj_atoms

   implicit none
   integer, intent(in) :: atom_of_func(:), atom_Z(:)

   integer              :: iatom, ifunc, f_count

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

end subroutine ljs_initialise

subroutine ljs_get_energy(energy)
   use LJ_switch_data, only: lj_atoms, mm_atoms

   implicit none
   LIODBLE, intent(out) :: energy



end subroutine ljs_get_energy


subroutine ljs_get_dEdQ()
   implicit none

end subroutine ljs_get_dEdQ

subroutine ljs_add_fock_terms(fock, energ, rho, S_matrix, n_of_func, &
                              atom_Z)
   use LJ_switch_data, only: n_lj_atoms, lj_atoms, ljs_initialised
   implicit none
   integer, intent(in)    :: n_of_func(:), atom_Z(:)
   LIODBLE, intent(in)    :: S_matrix(:,:), rho(:,:)
   LIODBLE, intent(inout) :: fock(:,:), energ

   integer :: ifunc, jfunc, iatom, f_idx
   LIODBLE :: atom_Q

   if (n_lj_atoms < 1) return

   if (.not. ljs_initialised) then
      call ljs_initialise(n_of_func, atom_Z)
      ljs_initialised = .true.
   endif

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

end module LJ_switch