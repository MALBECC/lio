!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine writegeom(Nsize,Vatnum,Vatpos,oUnit)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in)     :: Nsize
  integer,intent(in)     :: Vatnum(Nsize)
  real*8,intent(in)      :: Vatpos(3,Nsize)
  integer,intent(in)     :: oUnit

  integer :: kk


  write(Unit=oUnit,fmt=101) Nsize
  write(Unit=oUnit,fmt=102) ''
!  '  Atom Num        x             y             z'
!  xx[------]xx[----------]xx[----------]xx[----------]

  do kk=1,Nsize
    write(Unit=oUnit,fmt=103) Vatnum(kk), &
    Vatpos(1,kk)*(0.529177d0),            &
    Vatpos(2,kk)*(0.529177d0),            &
    Vatpos(3,kk)*(0.529177d0)
  enddo

101 format(6x,I4)
102 format(A)
103 format(2x,I4,3(2x,F12.8))

  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
