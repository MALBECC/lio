!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine matdcmp_svd( Matrix, Umat, Gmat, Vtrp )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   real*8 , intent(in)  :: Matrix(:,:)
   real*8 , intent(out) :: Umat(:,:)
   real*8 , intent(out) :: Gmat(:,:)
   real*8 , intent(out) :: Vtrp(:,:)

   logical              :: error_found
   integer              :: Msize1, Msize2
   real*8 , allocatable :: Xmat(:,:)
   real*8 , allocatable :: Svec(:)

   integer              :: lapack_LWORK
   real*8 , allocatable :: lapack_WORK(:)
   integer, allocatable :: lapack_IWORK(:)
   integer              :: lapack_INFO
!
!
!  Check sizes
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
      print*; stop
   endif
!
!
!  Obtains Lmat using lapack
!------------------------------------------------------------------------------!
   allocate( Xmat(Msize1,Msize2) )
   Xmat = Matrix

   allocate( Svec( min(Msize1,Msize2) ) )
   Svec = 0.0d0
!
!
!  QUERY
!------------------------------------------------------------------------------!
   lapack_LWORK = -1
   allocate( lapack_WORK(1) )
   allocate( lapack_IWORK( 8 * min(Msize1,Msize2) ) )
   call dgesdd( 'A', Msize1, Msize2, Xmat, Msize1, Svec, Umat, Msize1, Vtrp,   &
              & Msize2, lapack_WORK, lapack_LWORK, lapack_IWORK, lapack_INFO )

   if (lapack_INFO /= 0) then
      print*, 'ERROR INSIDE matdcmp_svd'
      print*, 'Problem in the call to dgesdd for query...'
      print*; stop
   endif
!
!
!  CALCULATION
!------------------------------------------------------------------------------!
   lapack_LWORK = NINT( lapack_WORK(1) )
   deallocate( lapack_WORK )
   allocate( lapack_WORK(lapack_LWORK) )
   call dgesdd( 'A', Msize1, Msize2, Xmat, Msize1, Svec, Umat, Msize1, Vtrp,   &
              & Msize2, lapack_WORK, lapack_LWORK, lapack_IWORK, lapack_INFO )

   if (lapack_INFO /= 0) then
      print*, 'ERROR INSIDE matdcmp_svd'
      print*, 'Problem in the call to dgesdd for calculations...'
      print*; stop
   endif

!
!
   deallocate( Xmat, Svec, lapack_WORK, lapack_IWORK )
end subroutine matdcmp_svd
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
