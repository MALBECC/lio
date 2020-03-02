!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_nucpos( nat, nucid0, nucpos, uid )
   implicit none
   integer         , intent(in) :: nat
   character(len=3), intent(in) :: nucid0(nat)
   LIODBLE    , intent(in) :: nucpos(nat,3)
   integer         , intent(in) :: uid

   character(len=*), parameter  :: fmtstr='(1X,A2,3(1X,F12.6))'
   integer :: kk

   write(unit=uid, fmt='(2X,I6)') nat
   write(unit=uid, fmt='(A)')
   do kk = 1, nat
       write(unit=uid, fmt=fmtstr) adjustl(nucid0(kk)), nucpos(kk,:)
   end do

end subroutine write_nucpos
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine translate_atomlist( natom, atom_idns, atom_idcs )
   implicit none
   integer         , intent(in)  :: natom
   integer         , intent(in)  :: atom_idns(natom)
   character(len=3), intent(out) :: atom_idcs(natom)
   integer                       :: kk

   do kk = 1, natom
       call atom_name( atom_idns(kk), atom_idcs(kk) )
   end do

end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
