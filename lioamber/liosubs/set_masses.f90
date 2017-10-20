!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine set_masses(Nsize,Vatnum,Vmass)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in)   :: Nsize
  integer,intent(in)   :: Vatnum(Nsize)
  real*8,intent(out)   :: Vmass(Nsize)
  integer              :: kk


  do kk=1,Nsize
    select case (Vatnum(kk))
      case (1);  Vmass(kk)=1.007940d0
      case (2);  Vmass(kk)=4.002602d0
      case (3);  Vmass(kk)=6.941000d0
      case (4);  Vmass(kk)=9.012182d0
      case (5);  Vmass(kk)=10.81100d0
      case (6);  Vmass(kk)=12.01070d0
      case (7);  Vmass(kk)=14.00670d0
      case (8);  Vmass(kk)=15.99940d0
      case (9);  Vmass(kk)=18.99840d0
      case (10); Vmass(kk)=20.17970d0
      case default
        print*,'ATOMIC NUMBER IS NOT SUPPORTED'
        stop
    endselect
!   Convert UMA to atomic mass:
    Vmass(kk)=Vmass(kk)*(1822.88857149)
  enddo


  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
