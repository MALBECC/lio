subroutine ObtainOsc(dip,E,O,N,Sref)
use excited_data, only: print_trdip
   implicit none

   integer, intent(in) :: N, Sref
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

   ! Print Transition Dipole Moment
   if ( print_trdip ) then
      open (unit=456,file="TransDipMom.dat")
      if ( Sref == 0 ) then
         write(456,*) "# Transtion Dipole Moments of GS->ES"
      else
         write(456,*) "# Transtion Dipole Moments of ES->ES"
      endif
      write(456,*) "# Transition    X        Y        Z    |Dip|^2 [a.u]"
      write(456,*) " "

      do ii=1,N
         temp = dip(ii,1)*dip(ii,1) + dip(ii,2)*dip(ii,2) + dip(ii,3)*dip(ii,3)
         write(456,"(I4,1X,I4,3X,F8.4,1X,F8.4,1X,F8.4,1X,F8.4)") Sref,Sref+ii,dip(ii,1),dip(ii,2),dip(ii,3), &
                                                         & temp
      enddo
   endif
end subroutine ObtainOsc
