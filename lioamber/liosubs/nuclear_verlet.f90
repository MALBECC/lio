!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine nuclear_verlet &
  (Npart,dt,mass,force,oldpos,nowpos,newpos,nowvel,Kenergy)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in)     :: Npart ! Degrees of Freedom
  real*8,intent(in)      :: dt  ! Time step

  real*8,intent(in)      :: mass(Npart)
  real*8,intent(in)      :: force(Npart,3)
  real*8,intent(in)      :: oldpos(Npart,3)
  real*8,intent(in)      :: nowpos(Npart,3)
  real*8,intent(out)     :: newpos(Npart,3)
  real*8,intent(out)     :: nowvel(Npart,3)
  real*8,intent(out)     :: Kenergy

  integer :: nn,kk
  real*8  :: dt2,dtsq


  dtsq=dt*dt
  dt2=2.0d0*dt
  Kenergy=0.0d0

  do kk=1,3
  do nn=1,Npart
    newpos(nn,kk)=2.0d0*nowpos(nn,kk)
    newpos(nn,kk)=newpos(nn,kk)-oldpos(nn,kk)
    newpos(nn,kk)=newpos(nn,kk)+(dtsq*force(nn,kk))/mass(nn)

    nowvel(nn,kk)=newpos(nn,kk)-oldpos(nn,kk)
    nowvel(nn,kk)=nowvel(nn,kk)/dt2
    Kenergy=Kenergy+(0.5d0)*mass(nn)*(nowvel(nn,kk)**2)
  enddo
  enddo

  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
