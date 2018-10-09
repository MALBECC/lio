!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine line_search( n_points, Energy, step_size, lambda )
!
!  If minimun value of E is E(n_points)
!     return lambda = step_size * n_points
!  If minimun value of E is E(1)
!     return lambda = 0.d0
!  else
!     returns an expected lambda that minimice E(lambda) using a
!     parabolic interpolation
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   integer         , intent(in)  :: n_points
   double precision, intent(in)  :: Energy(n_points)
   double precision, intent(in)  :: step_size
   double precision, intent(out) :: lambda

   integer          :: i
   integer          :: min_Energy_position
   double precision :: dE1, dE2, modif_fac

   if (n_points .le. 2) then
      write(*,*) "wrong n_points in lineal search, n_points need to be > 2"
      stop
   end if

!  find min value in Energy elements
   min_Energy_position=1
   do i=2, n_points
      if (Energy(i) .lt. Energy(min_Energy_position)) min_Energy_position=i
   end do

   if (min_Energy_position .eq. 1) then
     lambda=0.d0
     return
   elseif (min_Energy_position .eq. n_points) then
     lambda=step_size*dble(n_points)
     return
   end if

   dE2=abs(Energy(min_Energy_position) - Energy(min_Energy_position+1))
   dE1=abs(Energy(min_Energy_position) - Energy(min_Energy_position-1))

   modif_fac=step_size*(dE2-dE1)/(dE1+dE2)
   lambda=step_size*dble(min_Energy_position) - 0.5d0 * modif_fac

end subroutine line_search
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

