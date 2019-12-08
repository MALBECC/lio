subroutine vec_init(Vec,N,vecnum)
   implicit none

   integer, intent(in) :: N, vecnum
   double precision, intent(out) :: Vec(N,vecnum)

   integer :: ii

   Vec = 0.0D0
   do ii=1,vecnum
      Vec(ii,ii) = 1.0D0
   enddo
end subroutine vec_init
