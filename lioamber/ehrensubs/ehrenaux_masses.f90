!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrenaux_masses( Nsize, Vatnum, Vmass )
!------------------------------------------------------------------------------!
!
! DESCRIPTION
! Masses are in uma
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use constants_mod , only: massprot_elec
  implicit none
  integer,intent(in)     :: Nsize
  integer,intent(in)     :: Vatnum(Nsize)
  LIODBLE,intent(out)     :: Vmass(Nsize)

  integer :: kk
  
  
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
      case (11);  Vmass(kk)=22.9897d0
      case (12);  Vmass(kk)=24.305d0
      case (13);  Vmass(kk)=26.9815d0
      case (14);  Vmass(kk)=28.0855d0
      case (15);  Vmass(kk)=30.9738d0
      case (16);  Vmass(kk)=32.06500d0
      case (17);  Vmass(kk)=35.453d0
      case (18);  Vmass(kk)=39.948d0
      case (19);  Vmass(kk)=39.0983d0
      case (20);  Vmass(kk)=40.078d0
      case (21);  Vmass(kk)=44.9559d0
      case (22);  Vmass(kk)=47.867d0
      case (23);  Vmass(kk)=50.9415d0
      case (24);  Vmass(kk)=51.9961d0
      case (25);  Vmass(kk)=54.938d0
      case (26);  Vmass(kk)=55.84500d0
      case (27);  Vmass(kk)=58.9332d0
      case (28);  Vmass(kk)=59.6934d0
      case (29);  Vmass(kk)=63.546d0
      case (30);  Vmass(kk)=65.39d0
      case (31);  Vmass(kk)=69.723d0
      case (32);  Vmass(kk)=72.64d0
      case (33);  Vmass(kk)=74.9216d0
      case (34);  Vmass(kk)=78.96d0
      case (35);  Vmass(kk)=79.904d0
      case (36);  Vmass(kk)=83.8d0
      case (37);  Vmass(kk)=85.4678d0
      case (38);  Vmass(kk)=87.62d0
      case (39);  Vmass(kk)=88.9059d0
      case (40);  Vmass(kk)=91.224d0
      case (41);  Vmass(kk)=92.9064d0
      case (42);  Vmass(kk)=95.94d0
      case (43);  Vmass(kk)=98d0
      case (44);  Vmass(kk)=101.07d0
      case (45);  Vmass(kk)=102.9055d0
      case (46);  Vmass(kk)=106.42d0
      case (47);  Vmass(kk)=107.8682d0
      case (48);  Vmass(kk)=112.441d0
      case (49);  Vmass(kk)=114.818d0
      case (50);  Vmass(kk)=118.71d0
      case (51);  Vmass(kk)=121.76d0
      case (52);  Vmass(kk)=127.6d0
      case (53);  Vmass(kk)=126.9045d0
      case (54);  Vmass(kk)=131.293d0
      case default
        print*,'ATOMIC NUMBER ',Vatnum(kk),' IS NOT SUPPORTED'
        stop
    endselect
    Vmass(kk)=Vmass(kk)*massprot_elec
  enddo

end subroutine ehrenaux_masses
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
