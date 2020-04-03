subroutine print_sigma(sgm,nst,time)
use fstsh_data, only: current_state, tsh_file
   implicit none

   integer, intent(in) :: nst
   LIODBLE, intent(in) :: sgm(nst,nst)
   character(len=7), intent(in)   :: time

   integer :: ii, jj

   write(tsh_file,"(1X,A,A,A)") "Sigma[",trim(time),"]"
   do ii=1,nst
      do jj=1,nst
         write (tsh_file, "(F10.5)", ADVANCE="NO") sgm(ii,jj)
      enddo
      write(tsh_file,*) " "
   enddo
end subroutine print_sigma

subroutine print_Ener(Ene,istat,nstates)
use fstsh_data, only: tsh_file
   implicit none
   
   integer, intent(in) :: istat, nstates
   LIODBLE, intent(in) :: Ene(nstates)

   integer :: ii
   write(tsh_file,"(1X,A)",ADVANCE="NO") "Potential Energies= "
   do ii=1,nstates
      write (tsh_file, "(2X,F10.5)", ADVANCE="NO") Ene(ii)
   enddo
   write(tsh_file,*) " "
   write (tsh_file, "(1X,A,F10.5)") "Actual Potential Energy= ", Ene(istat)
end subroutine print_Ener
