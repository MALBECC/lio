!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine matdcmp_svd( Matrix, Umat, Gmat, Vtrp )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   LIODBLE , intent(in)  :: Matrix(:,:)
   LIODBLE , intent(out) :: Umat(:,:)
   LIODBLE , intent(out) :: Gmat(:,:)
   LIODBLE , intent(out) :: Vtrp(:,:)

   logical              :: error_found
   integer              :: nn
   integer              :: Msize1, Msize2
   LIODBLE , allocatable :: Xmat(:,:)
   LIODBLE , allocatable :: Gvec(:)

   integer              :: lapack_LWORK
   LIODBLE , allocatable :: lapack_WORK(:)
   integer, allocatable :: lapack_IWORK(:)
   integer              :: lapack_INFO
!
!
!  Checks and initialization
!------------------------------------------------------------------------------!
   Msize1 = size(Matrix,1)
   Msize2 = size(Matrix,2)
   error_found = .false.
   error_found = (error_found) .or. ( Msize1 /= size(Umat,1) )
   error_found = (error_found) .or. ( Msize1 /= size(Umat,2) )
   error_found = (error_found) .or. ( Msize1 /= size(Gmat,1) )
   error_found = (error_found) .or. ( Msize2 /= size(Gmat,2) )
   error_found = (error_found) .or. ( Msize2 /= size(Vtrp,1) )
   error_found = (error_found) .or. ( Msize2 /= size(Vtrp,2) )

   if (error_found) then
      print*, 'ERROR INSIDE matdcmp_svd'
      print*, 'Wrong sizes of input/output'
      print*,size(Matrix,1), size(Umat,1), size(Umat,2), size(Gmat,1)
      print*,size(Matrix,2), size(Vtrp,1), size(Vtrp,2), size(Gmat,2)
      print*; stop
   endif

   allocate( Xmat(Msize1,Msize2) )
   Xmat = Matrix

   allocate( Gvec( min(Msize1,Msize2) ) )
   Gvec(:) = 0.0d0
!
!
!  Query
!------------------------------------------------------------------------------!
   lapack_LWORK = -1
   allocate( lapack_WORK(1) )
   allocate( lapack_IWORK( 8 * min(Msize1,Msize2) ) )
   call dgesdd( 'A', Msize1, Msize2, Xmat, Msize1, Gvec, Umat, Msize1, Vtrp,   &
              & Msize2, lapack_WORK, lapack_LWORK, lapack_IWORK, lapack_INFO )

   if (lapack_INFO /= 0) then
      print*, 'ERROR INSIDE matdcmp_svd'
      print*, 'Problem in the call to dgesdd for query...'
      print*; stop
   endif
!
!
!  Calculation
!------------------------------------------------------------------------------!
   lapack_LWORK = INT( lapack_WORK(1) )
   deallocate( lapack_WORK )
   allocate( lapack_WORK(lapack_LWORK) )
   call dgesdd( 'A', Msize1, Msize2, Xmat, Msize1, Gvec, Umat, Msize1, Vtrp,   &
              & Msize2, lapack_WORK, lapack_LWORK, lapack_IWORK, lapack_INFO )

   if (lapack_INFO /= 0) then
      print*, 'ERROR INSIDE matdcmp_svd'
      print*, 'Problem in the call to dgesdd for calculations...'
      print*; stop
   endif
!
!
!  Copy results and exit
!------------------------------------------------------------------------------!
   Gmat(:,:) = 0.0d0
   do nn = 1, min(Msize1,Msize2)
      Gmat(nn,nn) = Gvec(nn)
   enddo

   deallocate( Xmat, Gvec, lapack_WORK, lapack_IWORK )
end subroutine matdcmp_svd
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
