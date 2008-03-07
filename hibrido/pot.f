c------------------------------------------------------------------
c subroutine for evaluating exchange correlation density and
c potential , depends of Iexch, index that says which potential
c will be used
c 1 X alpha
c 2 Gunnarson-Lundqvist
c 3 Vosko, Wilk and Nusair
c Exchange part is in all cases the same, for the time being LD
c Self interaction corrections are used in correlation part
c------------------------------------------------------------------
c
c     
      subroutine pot(Iexch,dens,ex,ec,v)
      implicit real*8 (a-h,o-z)
c
c data X alpha
      data const /-0.738558766382022447D0/
c data Gunnarson-Lundvquist
      data const2 /0.620350490899400087D0 /
c  data Vosko et al
      data A1,b1,c1,x0,Q,A16 /0.03109205D0,3.72744D0,12.9352D0,
     >  -0.10498D0,6.15199066246304849D0, 0.005182008333D0 /
      data A2,b2,c2,x02,Q2,A26 /0.015546025D0,7.06042D0,18.0578D0,
     > -0.32500D0,4.7309269D0,0.0025910042D0 /
c
        if (dens.eq.0.0D0) then
         ex=0.0D0
         ec=0.D0
         v = 0.D0 
         return
        endif
c       
        y=dens**0.333333333333333333D0
        e0=const*y
        v0=1.33333333333333D0*e0
c
      goto (10,20,30) Iexch
c
 10   continue
        ex=e0
        ec = 0.D0 
        v=v0
      return
c
 20   continue
        ex = e0
*
        rs=const2/y
        x1=rs/11.4D0
c
       if (x1.gt.1.D20) then
          ec=-0.0333D0*(0.5D0*x1-0.33333333333333D0)
          vc=0.0111D0*x1*0.500D0
       else
        t1=(1.D0+x1**3)
        t2=log(1.D0+1.D0/x1)
        t3=x1**2
        ec=-0.0333D0*(t1*t2-t3+0.5D0*x1-0.33333333333333D0)
        vc=0.0111D0*x1*(3.D0*t3*t2 - t1/(x1*(x1+1.D0))-2.D0*x1+0.5D0)
       endif
c self interaction correction, conceptually is not clear if it
c should be included or not
c       rs=rs*1.2599210585704356D0
c       x1=rs/15.9D0
c       t1=(1.D0+x1**3)
c       t2=log(1.D0+1.D0/x1)
c       t3=x1**2
c       ec2=-0.0203D0*(t1*t2-t3+0.5D0*x1-0.33333333333333D0)
c       vc2=0.0067666666D0*x1*(3.D0*t3*t2 - t1/(x1*(x1+1.D0))-
c    >      2.D0*x1+0.5D0)
c
c       e=e0+ec-ec2
c       v=v0+ec+vc-vc2
*        e=e0+ec
        v=v0+ec+vc
        return
c
 30   continue
c
        ex = e0
*
        rs=const2/y
        x1=sqrt(rs)
        Xx=rs+b1*x1+c1
        Xxo=x0**2+b1*x0+c1
        t1=2.D0*x1+b1
        t2=log(Xx)
        t3=atan(Q/t1)
        t4=b1*x0/Xxo
c
        ec=A1*(2.D0*log(x1)- t2+ 2.D0*b1/Q*t3 -t4*(2.D0*log(x1-x0)-
     >     t2+ 2.D0*(b1+2.D0*x0)/Q*t3))
c 
        t5=(b1*x1+2.D0*c1)/x1
        t6=x0/Xxo
        vc=ec-A16*x1*(t5/Xx-4.D0*b1/(t1**2+Q**2)*(1.D0-t6*(b1-2.D0*x0))
     >     -t4*(2.D0/(x1-x0)-t1/Xx))
c
c self interaction correction
c       rs=rs*1.2599210585704356D0
c       x1=sqrt(rs)
c       Xx=rs+b2*x1+c2
c       Xxo=x02**2+b2*x02+c2
c       t1=2.D0*x1+b2
c       t2=log(Xx)
c       t3=atan(Q2/t1)
c       t4=b2*x02/Xxo
c
c       e2=A2*(2.D0*log(x1)- t2+ 2.D0*b2/Q2*t3 -t4*(2.D0*log(x1-x02)-
c    >     t2+ 2.D0*(b2+2.D0*x02)/Q2*t3))
c
c       t5=b2*x1**2+2.D0*c2*x1
c       t6=x02/Xxo
c      v2=e2-A26*x1*(t5/Xx-4.D0*b2/(t1**2+Q2**2)*(1.D0-t6*(b2-2.D0*x02))
c    >    -t4*(2.D0/(x1-x02)-t1/Xx))
c
c       e=e0+ec-e2
c       v=v0+vc-v2
c
*        e=e0+ec
        v=v0+vc
        return
c
      end
