subroutine ObtainOsc(dip,E,O,N)
   implicit none

   integer, intent(in) :: N
   LIODBLE, intent(in)  :: dip(N,3), E(N)
   LIODBLE, intent(out) :: O(N)

   integer :: ii, jj
   LIODBLE :: dostres, temp

   dostres = 2.0D0 / 3.0D0
   temp = 0.0D0

   do ii=1,N
   do jj=1,3
      temp = temp + dip(ii,jj) * dip(ii,jj) * E(ii) * dostres
   enddo
   O(ii) = temp; temp = 0.0D0
   enddo
end subroutine ObtainOsc
