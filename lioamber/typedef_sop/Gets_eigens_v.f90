!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Gets_eigens_v( this, maxval_ld, eigenvals, eigenvecs )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   implicit none
   class(sop), intent(in)            :: this
   real*8    , intent(in)            :: maxval_ld
   real*8    , intent(out)           :: eigenvals(:)
   real*8    , intent(out), optional :: eigenvecs(:,:)

   real*8    , allocatable           :: eigenvals_m(:,:)
   logical                           :: error_found
   integer                           :: nn


!  Checks and preps
!------------------------------------------------------------------------------!
   if ( this%Nbasis <= 0 ) then
      print*, "ERROR INSIDE Gets_eigens_v: overlap matrix was never set"
      print*, "ABORTING RUN"; stop
   end if

   error_found = .false.
   error_found = (error_found).or.( this%Nbasis /= size(eigenvals) )


!  Calculations
!------------------------------------------------------------------------------!
   allocate( eigenvals_m( this%Nbasis, this%Nbasis ) )
   if ( present(eigenvecs) ) then
      call Gets_eigens_m( this, maxval_ld, eigenvals_m, eigenvecs )
   else
      call Gets_eigens_m( this, maxval_ld, eigenvals_m )
   end if
   do nn = 1, this%Nbasis
      eigenvals(nn) = eigenvals_m(nn,nn)
   end do
   deallocate( eigenvals_m )


end subroutine Gets_eigens_v
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
