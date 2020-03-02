!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine nuclear_verlet (Natom,dt,masses,forces,oldpos,nowpos,&
                             newpos,nowvel,Kenergy)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in)  :: Natom ! Number of atoms
  LIODBLE,intent(in)   :: dt    ! Time step

  LIODBLE,intent(in)   :: masses(Natom)
  LIODBLE,intent(in)   :: forces(Natom,3)
  LIODBLE,intent(in)   :: oldpos(Natom,3)
  LIODBLE,intent(in)   :: nowpos(Natom,3)

  LIODBLE,intent(out)  :: newpos(Natom,3)
  LIODBLE,intent(out)  :: nowvel(Natom,3)
  LIODBLE,intent(out)  :: Kenergy

  integer :: nn,kk
  LIODBLE  :: dt2,dtsq

!--------------------------------------------------------------------!

  dtsq=dt*dt
  dt2=2.0d0*dt
  Kenergy=0.0d0

  do kk=1,3
  do nn=1,Natom
    newpos(nn,kk)=(2.0d0)*nowpos(nn,kk)
    newpos(nn,kk)=newpos(nn,kk)-oldpos(nn,kk)
    newpos(nn,kk)=newpos(nn,kk)+(dtsq*forces(nn,kk))/masses(nn)

    nowvel(nn,kk)=newpos(nn,kk)-oldpos(nn,kk)
    nowvel(nn,kk)=nowvel(nn,kk)/dt2
    Kenergy=Kenergy+(0.5d0)*masses(nn)*(nowvel(nn,kk)**2)
  enddo
  enddo

  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
