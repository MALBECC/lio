!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrenaux_calckyn( Npart, mass, vels, Kenergy )
!------------------------------------------------------------------------------!
!
! DESCRIPTION
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in)     :: Npart
  LIODBLE,intent(in)      :: mass(Npart)
  LIODBLE,intent(in)      :: vels(3,Npart)
  LIODBLE,intent(out)     :: Kenergy
  integer :: nn,kk

  Kenergy=0.0d0
  do nn=1,Npart
  do kk=1,3
    Kenergy=Kenergy+(0.5d0)*mass(nn)*(vels(kk,nn)**2)
  enddo
  enddo
  Kenergy=Kenergy*(0.0005485799d0)

end subroutine ehrenaux_calckyn
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
