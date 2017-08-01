!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
! Save the forces and density matrix of the last position for restart purposes
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rstsave( rstunit, Natom, forces, Nbasis, densA, densB )
   implicit none
   integer,    intent(in) :: rstunit
   integer,    intent(in) :: Natom
   real*8,     intent(in) :: forces( 3, Natom )
   integer,    intent(in) :: Nbasis
   complex*16, intent(in) :: densA( Nbasis, Nbasis )
   complex*16, intent(in) :: densB( Nbasis, Nbasis )

   integer :: ii, jj

   do jj=1,Natom
   do ii=1,3
      write( unit=rstunit, fmt=100 ) ii, jj, forces(ii, jj)
   enddo
   enddo

   do jj=1,Nbasis
   do ii=1,Nbasis
      write( unit=rstunit, fmt=100 ) ii, jj, densA(ii, jj), densB(ii, jj)
   enddo
   enddo

100 format(2x,I3,2x,I3,4(2x,ES20.12))
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rstload( rstunit, Natom, forces, Nbasis, densA, densB )
   implicit none
   integer,    intent(in)  :: rstunit
   integer,    intent(in)  :: Natom
   real*8,     intent(out) :: forces( 3, Natom )
   integer,    intent(in)  :: Nbasis
   complex*16, intent(out) :: densA( Nbasis, Nbasis )
   complex*16, intent(out) :: densB( Nbasis, Nbasis )

   integer :: ii, jj, isc, jsc

   do jj=1,Natom
   do ii=1,3
      read( unit=rstunit, fmt=100 ) isc, jsc, forces(ii, jj)
   enddo
   enddo

   do jj=1,Nbasis
   do ii=1,Nbasis
      read( unit=rstunit, fmt=100 ) isc, jsc, densA(ii, jj), densB(ii,jj)
   enddo
   enddo

100 format(2x,I3,2x,I3,4(2x,ES20.12))
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
