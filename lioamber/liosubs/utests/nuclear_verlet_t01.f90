!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  include "../nuclear_verlet.f90"
  program nuclear_verlet_t01
  implicit none
  real*8  :: oldpstn(2,3), nowpstn(2,3), newpstn(2,3)
  real*8  :: nowvels(2,3), nowfrce(2,3)
  real*8  :: Kenergy,Uenergy,fconst,mass(2),dt
  real*8  :: posvec(3), dist, dist_eq
  integer :: nn,kk,nsteps


! Setups
!------------------------------------------------------------------------------!
  nsteps=10200
  dt=0.1d0
  dist_eq=0.5d0
  fconst=1.0d0
  mass(1)=1.0d0
  mass(2)=1.0d0
  nowpstn(1,1)=-0.2d0
  nowpstn(2,1)= 0.2d0
  nowpstn(1,2)=-0.2d0
  nowpstn(2,2)= 0.2d0
  nowpstn(1,3)=-0.2d0
  nowpstn(2,3)= 0.2d0

! Initialization
!------------------------------------------------------------------------------!
  oldpstn=nowpstn
  oldpstn(2,1)= 0.199d0
  open(unit=101,file='Output-nuclear_verlet_t01_E.dat')
  open(unit=102,file='Output-nuclear_verlet_t01_Q.xyz')
  write(unit=102,fmt=202) '    2'
  write(unit=102,fmt=202) ''
  write(unit=102,fmt=202) '1',nowpstn(1,1), nowpstn(1,2), nowpstn(1,3)
  write(unit=102,fmt=202) '1',nowpstn(2,1), nowpstn(2,2), nowpstn(2,3)


! Dinamic Steps
!------------------------------------------------------------------------------!
  do nn=1,nsteps

    dist=0.0d0
    do kk=1,3
      posvec(kk)=nowpstn(1,kk)-nowpstn(2,kk)
      dist=dist+posvec(kk)**2
    enddo
    dist=sqrt(dist)
    Uenergy=(0.5)*(dist-dist_eq)**2

    do kk=1,3
      nowfrce(1,kk)=(1.0d0)-(dist_eq/dist)
      nowfrce(1,kk)=nowfrce(1,kk)*posvec(kk)
      nowfrce(1,kk)=nowfrce(1,kk)*(-1.0d0)
      nowfrce(1,kk)=nowfrce(1,kk)*(fconst)

      nowfrce(2,kk)=nowfrce(1,kk)*(-1.0d0)
    enddo

    call nuclear_verlet &
      (2,dt,mass,nowfrce,oldpstn,nowpstn,newpstn,nowvels,Kenergy)
    oldpstn=nowpstn
    nowpstn=newpstn

    write(unit=101,fmt=201) nn,Uenergy+Kenergy,Uenergy,Kenergy
    write(unit=102,fmt=202) '    2'
    write(unit=102,fmt=202) ''
    write(unit=102,fmt=202) '1',nowpstn(1,1), nowpstn(1,2), nowpstn(1,3)
    write(unit=102,fmt=202) '1',nowpstn(2,1), nowpstn(2,2), nowpstn(2,3)
  enddo

!------------------------------------------------------------------------------!
  201 format(I6,3(2x,f10.6))
  202 format(A,3(2x,f10.6))
  return; end program
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
