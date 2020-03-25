subroutine print_sigma(sgm,nst)
use fstsh_data, only: current_state, tsh_file
   implicit none

   integer, intent(in) :: nst
   LIODBLE, intent(in) :: sgm(nst,nst)

   integer :: ii, jj

   write(tsh_file,"(1X,A)") "Sigma"
   do ii=1,nst
      do jj=1,nst
         write (tsh_file, "(F10.5)", ADVANCE="NO") sgm(ii,jj)
      enddo
      write(tsh_file,*) " "
   enddo
end subroutine print_sigma
