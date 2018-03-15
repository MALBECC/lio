!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Gets_eigens_m( this, maxval_ld, eigenvals, eigenvecs )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   implicit none
   class(sop), intent(in)            :: this
   real*8    , intent(in)            :: maxval_ld
   real*8    , intent(out)           :: eigenvals(:,:)
   real*8    , intent(out), optional :: eigenvecs(:,:)

   logical                           :: error_found
   integer                           :: nn


!  Checks and preps
!------------------------------------------------------------------------------!
   if ( this%Nbasis <= 0 ) then
      print*, "ERROR INSIDE Gets_eigens_m: overlap matrix was never set"
      print*, "ABORTING RUN"; stop
   end if

   error_found = .false.
   error_found = (error_found).or.( this%Nbasis /= size(eigenvals,1) )
   error_found = (error_found).or.( this%Nbasis /= size(eigenvals,2) )
   if ( present(eigenvecs) ) then
      error_found = (error_found).or.( this%Nbasis /= size(eigenvecs,1) )
      error_found = (error_found).or.( this%Nbasis /= size(eigenvecs,2) )
   end if


!  Calculations
!------------------------------------------------------------------------------!
   if ( present(eigenvecs) ) eigenvecs = this%Umat

   eigenvals = this%Gmat
   call this%Drop_ldvals( maxval_ld, eigenvals )
   do nn = 1, this%Nbasis
      eigenvals(nn,nn) = eigenvals(nn,nn) * eigenvals(nn,nn)
   end do

end subroutine Gets_eigens_m
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
