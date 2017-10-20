!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine write_energy(Time,Epot,Ekin,Etot,oUnit)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  real*8,intent(in)  :: Time
  real*8,intent(in)  :: Epot
  real*8,intent(in)  :: Ekin
  real*8,intent(in)  :: Etot
  integer,intent(in) :: oUnit

  real*8 :: Time_ps


  Time_ps=Time*0.000024189
  if (Time.lt.0) then
     write(Unit=oUnit,fmt=101) &
    '   Time (ps)     Potential (au)   Kinetic (au)     Total Energy (au)'
!    xx[---------]xxxx[------------]xx[------------]xxxxxx[----------]
  else
     write(Unit=oUnit,fmt=102) Time_ps,Epot,Ekin,Etot
  endif


101 format(A)
102 format(2x,F11.7,2x,2(2x,F14.8),6x,F14.8)
  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
