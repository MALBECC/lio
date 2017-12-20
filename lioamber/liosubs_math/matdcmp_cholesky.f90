!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine matdcmp_cholesky( Matrix, Lmat )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   real*8 , intent(in)  :: Matrix(:,:)
   real*8 , intent(out) :: Lmat(:,:)

   integer :: Msize, ii, jj
   logical :: error_found
   integer :: lapack_INFO
!
!
!  Check sizes
!------------------------------------------------------------------------------!
   Msize = size(Matrix,1)
   error_found = .false.
   error_found = (error_found) .or. ( Msize /= size(Matrix,2) )
   error_found = (error_found) .or. ( Msize /= size(Lmat,1) )
   error_found = (error_found) .or. ( Msize /= size(Lmat,2) )

   if (error_found) then
      print*, 'ERROR INSIDE matdcmp_cholesky'
      print*, 'Wrong sizes of input/output'
      print*; stop
   endif
!
!
!  Obtains Lmat using lapack
!------------------------------------------------------------------------------!
   Lmat = Matrix
   call dpotrf( 'L', Msize, Lmat, Msize, lapack_INFO )
   if (lapack_INFO /= 0) then
      print*, 'ERROR INSIDE matdcmp_cholesky'
      print*, 'Problem in the call to dpotrf...'
      print*; stop
   endif

   do ii = 1, Msize-1
   do jj = ii+1, Msize
      Lmat(ii,jj) = 0.0d0
   enddo
   enddo
!
!
end subroutine matdcmp_cholesky
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
