!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% RESTART_FOCK.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! * read_fock_restart_c   (reads Fock restart for closed shell)                !
! * read_fock_restart_o   (reads Fock restart for open shell)                  !
! * write_fock_restart_c  (writes Fock restart for closed shell)               !
! * write_fock_restart_o  (writes Fock restart for open shell)                 !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

subroutine write_fock_restart_c(fock, M, n_occ, ndx, UID)
   ! fock  : Coefficient matrix.
   ! M     : Total number of basis functions.
   ! n_occ : Total number of occupied orbitals.
   ! ndx   : Basis functions index.
   ! UID   : Output file unit.
   implicit none
   integer, intent(in) :: M, n_occ, UID, ndx(M)
   real*8 , intent(in) :: fock(M,4*M)

   integer :: icount, jcount

   do icount = 1, M
   do jcount = 1, M
      fock(ndx(icount), M+jcount) = fock(icount, 2*M+jcount)      
   enddo
   enddo

   rewind UID
   do icount = 1, M
   do jcount = 1, n_occ
      write(UID, 400) fock(icount, M+jcount)
   enddo
   enddo   

   return

400 format(4(E14.7E2, 2x))
end subroutine write_fock_restart_c


subroutine write_fock_restart_o(fock_a, fock_b, M, n_occ_a, n_occ_b, ndx, UID)
   ! fock_a  : Alpha coefficient matrix.
   ! fock_b  : Beta coefficient matrix.
   ! M       : Total number of basis functions.
   ! n_occ_a : Total number of occupied alpha orbitals.
   ! n_occ_a : Total number of occupied beta orbitals.
   ! ndx     : Basis functions index.
   ! UID     : Output file unit.
   implicit none
   integer, intent(in) :: M, n_occ_a, n_occ_b, UID, ndx(M)
   real*8 , intent(in) :: fock_a(M,4*M), fock_b(M, 4*M)

   integer :: icount, jcount

   do icount=1, M
   do jcount=1, M
      fock_a(ndx(icount), M+jcount) = fock_a(icount, 2*M+jcount)
      fock_b(ndx(icount), M+jcount) = fock_b(icount, 2*M+jcount)
   enddo
   enddo

   rewind UID
   do icount = 1, M
   do jcount = 1, n_occ_a
      write(UID, 400) fock_a(icount, M+jcount)
   enddo
   enddo

   do icount = 1, M
   do jcount = 1, n_occ_b
      write(UID, 400) fock_b(icount, M+jcount)
   enddo
   enddo

   return

400 format(4(E14.7E2, 2x))
end subroutine write_fock_restart_o


