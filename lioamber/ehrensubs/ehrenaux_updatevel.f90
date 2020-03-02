!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrenaux_updatevel( natoms, masses, forces, nucvel, dt )

   implicit none
   integer, intent(in)    :: natoms
   LIODBLE,  intent(in)    :: masses(natoms)
   LIODBLE,  intent(in)    :: forces(3, natoms)
   LIODBLE,  intent(inout) :: nucvel(3, natoms)
   LIODBLE,  intent(in)    :: dt

   integer                :: nn, kk

   do nn = 1, natoms
   do kk = 1, 3
      nucvel(kk,nn) = nucvel(kk,nn) + dt * forces(kk,nn) / masses(nn)
   enddo
   enddo

end subroutine ehrenaux_updatevel
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
