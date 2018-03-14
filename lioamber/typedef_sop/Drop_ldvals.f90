!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Drop_ldvals( maxval_ld, eigen_valsq, eigen_invsq )
!------------------------------------------------------------------------------!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   real*8    , intent(in)              :: maxval_ld
   real*8    , intent(inout)           :: eigen_valsq(:,:)
   real*8    , intent(inout), optional :: eigen_invsq(:,:)

   integer                             :: Msize, nn
   logical                             :: error_found

   Msize = size( eigen_valsq, 1 )
   error_found = .false.
   error_found = (error_found) .or. ( Msize /= size(eigen_valsq,2) )
   if ( present(eigen_invsq) ) then
      error_found = (error_found) .or. ( Msize /= size(eigen_invsq,1) )
      error_found = (error_found) .or. ( Msize /= size(eigen_invsq,2) )
   end if

   do nn = 1, Msize
      if ( eigen_valsq(nn,nn)*eigen_valsq(nn,nn) < maxval_ld ) then
         eigen_valsq(nn,nn) = 0.0d0
         if ( present(eigen_invsq) ) eigen_invsq(nn,nn) = 0.0d0
      end if
   end do

end subroutine Drop_ldvals
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
