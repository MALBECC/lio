
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
         write(*,print_fmt) "Atom: ", iatom, "Index: ", lj_atoms(iatom)%idx

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