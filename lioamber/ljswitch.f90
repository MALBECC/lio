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
   end type lj_atom

   integer :: n_lj_atoms
   type(lj_atom), allocatable :: lj_atoms(:)

end module LJ_switch_data

module LJ_switch
   implicit none
   private
   public :: lj_input_read

contains

subroutine lj_input_read(input_UID)
   use LJ_switch_data, only: n_lj_atoms, lj_atoms
   implicit none
   integer, intent(in) :: input_UID
    
   character(len=10) :: buffer
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

   if (allocated(lj_atoms)) deallocate(lj_atoms(n_lj_atoms))

   ! Starts reading LJ switch data.
   write(*,'(A)') " AtomID |  Q1  |  Q2   |  Sigma1  | Sigma2 | Epsilon1 "&
                 &"| Epsilon2 "
   do iatom = 1, n_lj_atoms
      read(input_UID,*) lj_atoms(iatom)%idx, lj_atoms(iatom)%q1, &
                        lj_atoms(iatom)%q2 , lj_atoms(iatom)%s1, &
                        lj_atoms(iatom)%s2 , lj_atoms(iatom)%e1, &
                        lj_atoms(iatom)%e2
   enddo
   
 end subroutine lj_input_read



end module LJ_switch