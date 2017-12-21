!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine purge_zeros_v( cutval, vector )

   implicit none
   real*8 , intent(in)    :: cutval
   real*8 , intent(inout) :: vector(:)
   integer                :: ii

   do ii = 1, size(vector,1)
      if ( abs(vector(ii)) < cutval ) vector(ii) = 0.0d0
   end do

end subroutine purge_zeros_v
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine purge_zeros_m( cutval, matrix )

   implicit none
   real*8 , intent(in)    :: cutval
   real*8 , intent(inout) :: matrix(:,:)
   integer                :: ii, jj

   do jj = 1, size(matrix,1)
   do ii = 1, size(matrix,2)
      if ( abs(matrix(ii,jj)) < cutval ) matrix(ii,jj) = 0.0d0
   end do
   end do

end subroutine purge_zeros_m
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
