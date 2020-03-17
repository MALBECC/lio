!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine gaussian_shaper( do_grow, do_fall, xpos, center, sigma1, yval )
   implicit none
   logical, intent(in)    :: do_grow, do_fall
   LIODBLE , intent(in)    :: xpos
   LIODBLE , intent(in)    :: center
   LIODBLE , intent(in)    :: sigma1
   LIODBLE , intent(inout) :: yval

   logical :: is_growing, is_falling
   LIODBLE  :: gauss_exp

   is_growing = (do_grow) .and. ( xpos < center )
   is_falling = (do_fall) .and. ( xpos > center )

   if ( (is_growing) .or. (is_falling) ) then
      gauss_exp = ( xpos - center ) / (sigma1)
      gauss_exp = (-1.0d0) * (gauss_exp)**2
      yval = yval * exp( gauss_exp )
   endif

end subroutine gaussian_shaper
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

