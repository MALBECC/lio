!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine write_geom(Nsize,Vatnum,Vatpos,oUnit)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in)     :: Nsize
  integer,intent(in)     :: Vatnum(Nsize)
  real*8,intent(in)      :: Vatpos(Nsize,3)
  integer,intent(in)     :: oUnit

  integer :: kk


  write(Unit=oUnit,fmt=101) Nsize
  write(Unit=oUnit,fmt=100) &
  '  Atom Num        x             y             z'
!  xxxx[--]xxxx[----------]xx[----------]xx[----------]

  do kk=1,Nsize
    write(Unit=oUnit,fmt=103) Vatnum(kk),                &
                              Vatpos(kk,1)*(0.529177d0), &
                              Vatpos(kk,2)*(0.529177d0), &
                              Vatpos(kk,3)*(0.529177d0)
  enddo

100 format(A)
101 format(6x,I4)
103 format(4x,I4,2x,3(2x,F12.8))

  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
