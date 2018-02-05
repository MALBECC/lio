!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine standard_coefs( coef_mat )
!
! dens_mat:   matrix containing output density matrix.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   real*8 , intent(inout) :: coef_mat(:, :)
   integer                :: ii, jj

   do jj = 1, size(coef_mat, 2)
      if ( coef_mat(1,jj) < 0.0d0 ) then
         do ii = 1, size(coef_mat, 1)
            coef_mat(ii,jj) = (-1.0d0) * coef_mat(ii,jj)
         end do
      end if
   end do

end subroutine standard_coefs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
