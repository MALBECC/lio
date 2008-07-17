      implicit real*16 (a-h,o-z)
      pi=acos(-1.0D0)
      cf=3.0D0/10.0D0*(3.0D0*pi**2)**(2.0D0/3.0D0)
      write(*,100) cf
 100  format(F11.8)
      end
