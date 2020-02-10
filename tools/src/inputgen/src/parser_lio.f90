!##############################################################################!
module parser_lio
   implicit none
contains


!##############################################################################!

subroutine writelio_zyx( funit, natoms, numids, coords )
   implicit none
   integer, intent(in) :: funit
   integer, intent(in) :: natoms
   integer, intent(in) :: numids(natoms)
   real*8,  intent(in) :: coords(3, natoms)
   integer             :: ii
   
   do ii = 1, natoms
      write( unit=funit, fmt="(I5,4F16.8)" ) numids(ii), coords(:, ii)
   end do
   
end subroutine writelio_zyx

!------------------------------------------------------------------------------!

subroutine writelio_rst( funit, mdim, densmat )
   implicit none
   integer, intent(in) :: funit
   integer, intent(in) :: mdim
   real*8,  intent(in) :: densmat(mdim, mdim)
   real*8, allocatable :: temp(:)
   integer             :: ii

   allocate( temp(mdim) )

   do ii = 1, mdim
      temp(:) = densmat(ii,:) * 2.0d0
      temp(ii) = densmat(ii,ii)
!      write( unit=funit, fmt='(4(E14.7E2, 2x))' ) densmat(ii,:)
      write( unit=funit, fmt=* ) temp
   end do

end subroutine writelio_rst

!##############################################################################!
end module parser_lio
!##############################################################################!
