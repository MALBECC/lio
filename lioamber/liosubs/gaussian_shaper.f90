!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine gaussian_shaper( do_grow, do_fall, xpos, center, sigma1, yval )
   implicit none
   logical, intent(in)    :: do_grow, do_fall
   real*8 , intent(in)    :: xpos
   real*8 , intent(in)    :: center
   real*8 , intent(in)    :: sigma1
   real*8 , intent(inout) :: yval

   logical :: is_growing, is_falling
   real*8  :: gauss_exp

   is_growing = (do_grow) .and. ( xpos < center )
   is_falling = (do_fall) .and. ( xpos > center )

   if ( (is_growing) .or. (is_falling) ) then
      gauss_exp = ( xpos - center ) / (sigma1)
      gauss_exp = (-1.0d0) * (gauss_exp)**2
      yval = yval * exp( gauss_exp )
   endif

end subroutine gaussian_shaper
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

