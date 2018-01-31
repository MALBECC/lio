!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine messup_densmat( dens_mat )
!
! dens_mat:   matrix containing output density matrix.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   real*8 , intent(inout) :: dens_mat(:, :)
   integer                :: ii, jj

   do jj = 1, size(dens_mat, 2)
   do ii = 1, size(dens_mat, 1)
      if ( ii /= jj ) dens_mat(ii, jj) = dens_mat(ii, jj) * 2.0d0
   enddo
   enddo

end subroutine messup_densmat
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
