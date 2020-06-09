#include "datatypes/datatypes.fh"
module LJ_switch_data
   implicit none

   type lj_atom
      integer :: idx
      LIODBLE :: q1
      LIODBLE :: q2
      LIODBLE :: s1
      LIODBLE :: s2
      LIODBLE :: e1
      LIODBLE :: e2
      LIODBLE :: eps
      LIODBLE :: sig

      contains
         procedure, pass :: set_eps_sig
   end type lj_atom

   integer :: n_lj_atoms
   type(lj_atom), allocatable :: lj_atoms(:)

contains

   subroutine set_eps_sig(this, chrg)
      implicit none
      LIODBLE       , intent(in)    :: chrg
      class(lj_atom), intent(inout) :: this

      LIODBLE :: exp_term

      exp_term = exp( -10.0D0 * ( chrg - 0.5D0 * (this%q2 + this%q1) ) )
      exp_term = 1.0D0 / (1.0D0 + exp_term)

      this%sig = this%s1 + exp_term * (this%s2 - this%s1)
      this%eps = this%e1 + exp_term * (this%e2 - this%e1)
   end subroutine set_eps_sig

end module LJ_switch_data

module LJ_switch
   implicit none
   private
   public :: doing_ljs
   public :: ljs_input_read
   public :: ljs_calc_params
   public :: ljs_get_params

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
   integer, intent(in) :: input_UID
   integer, intent(in) :: verbose_lvl
    
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
   
end subroutine ljs_input_read

subroutine ljs_calc_params(atom_crg)
   use LJ_switch_data, only: n_lj_atoms, lj_atoms

   implicit none
   LIODBLE, intent(in) :: atom_crg(:)
   integer :: iatom
   
   if (n_lj_atoms < 1) return

   do iatom = 1, n_lj_atoms
      call lj_atoms(iatom)%set_eps_sig( atom_crg(lj_atoms(iatom)%idx) )
   enddo

end subroutine ljs_calc_params

subroutine ljs_get_params(sig, eps, dsig, deps)
   use LJ_switch_data, only: n_lj_atoms, lj_atoms

   implicit none
   LIODBLE, intent(inout) :: sig(n_lj_atoms), eps(n_lj_atoms)
   LIODBLE, intent(inout) :: dsig(n_lj_atoms), deps(n_lj_atoms)
   integer :: iatom

   if (n_lj_atoms < 1) return
   do iatom = 1, n_lj_atoms
      sig(iatom) = lj_atoms(iatom)%sig
      eps(iatom) = lj_atoms(iatom)%eps
   enddo
end subroutine ljs_get_params

end module LJ_switch