
c------------------------------------------------------------------
c subroutine for evaluating exchange correlation density and
c potential , for non local density functionals
c Iexch, index should be greater than 4
c 4 : Exchange given by Perdew : Phys. Rev B 33 8800 (1986)
c     Correlation :     Perdew : Phys. Rev B 33 8822 (1986)
c
c 5 : Exchange given by Becke  : Phys. Rev A 38 3098 (1988)
c     Correlation :     Perdew : Phys. Rev B 33 8822 (1986) 
c 6 : only  exchange given by Becke
c to be used in concert with dnsg subroutine, the one that gives
c gradient and second derivatives of the electronic density
c 19-1-93
c------------------------------------------------------------------
c
c     
      subroutine potg(Iexch,dens,dx,dy,dz,dxx,dyy,dzz,dxy,dyz,dxz,
     >                ex,ec,v)
      implicit real*8 (a-h,o-z)
      external asinh
c
c data X alpha
      data const /-0.738558766382022447D0/
c data Gunnarson-Lundvquist
      data const2 /0.620350490899400087D0 /
c  data Vosko et al
      data A1,b1,c1,x0,Q,A16 /0.03109205D0,3.72744D0,12.9352D0,
     >  -0.10498D0,6.15199066246304849D0, 0.0051820083D0 /
      data A2,b2,c2,x02,Q2,A26 /0.015546025D0,7.06042D0,18.0578D0,
     >  -0.32500D0,4.7309269D0,0.0025910042D0 /
*  data Becke
      data beta /0.0042D0/
c
c non local Perdew functional
c 
      data alf,bet,gam,del /0.023266D0,7.389D-06,8.723D0,0.472D0/
c
        if (dens.eq.0) then
         v=0.0D0
         ex=0.0D0
         ec = 0.D0
         return
        endif
c
        y=dens**0.333333333333333333D0
c
c
c-- Non local --------------------------- 
c
c
c Perdew (Phys. Rev B, 33 8800)
c results checked for energies against paper
       ckF=3.0936677D0*y
       DGRAD2=Dx**2+Dy**2+Dz**2
       DGRAD=sqrt(DGRAD2)
       u0=(Dx**2*Dxx+ 2.D0*Dx*Dy*Dxy+ 2.D0*Dy*Dz*Dyz+ 2.D0*Dx*Dz*Dxz+
     >   Dy**2.0D0*Dyy+ Dz**2.0D0*Dzz)/DGRAD
       D0=Dxx+Dyy+Dzz 
*
        if(iexch.ge.5) goto 451
*
*     VWN's exchange     
*
       dens2=dens**2
       s=DGRAD/(2.D0*ckF*dens)
c
c   
       fx=1.D0/15.D0
       s2=s**2
       s3=s**3
       s4=s**4
       g0=1.D0+1.296D0*s2+14.D0*s4+0.2D0*s**6
       F=g0**fx
       e=const*F*y
       ex = e 
c
c
       t=D0/(dens*4.D0*ckF**2)
       u=u0/((2.D0*ckF)**3*dens2)
c
       g2=2.592D0*s+56.D0*s3+1.2D0*s**5
       g3=2.592D0+56.D0*s2+1.2D0*s4
       g4=112.D0*s+4.8D0*s3
       dF=fx*F/g0*g2
       dsF=fx*F/g0*(-14.D0*fx*g3*g2/g0+g4)
       v=const*(1.33333333333D0*F-t/s*dF-(u-1.3333333333D0*s3)*dsF)*y
       goto 452
c In some cases, instabilities occur, and an ad-hoc solution is to use
c this expression for the potential, in those cases the exchange 
c potential of Becke should be used 
c      v=const*1.33333333333D0*F*y
c -----
*********
*        Becke's exchange
*
****
*        energy
 451   e = 0.0D0
       v = 0.0D0
       e0 = const*y
       y2 = dens/2.0D0
       r13 = y2**0.333333333333333333D0
       r43 = y2*r13
       Xs = DGRAD/(2.0D0*r43)
       siper = asinh(Xs)
       DN = 1.0D0 + 6.0D0*beta*Xs*siper
       ect = - 2.0D0*beta*r43*Xs*Xs/(DN*dens)
       e = e0 + ect 
       ex = e
*
*       return
*        potential
       v0 = 1.33333333333333D0*e0 
       Fb = 1.0D0/DN
       XA1 = Xs/dsqrt(1+Xs*Xs)
       DN1 = 1.0D0 + Fb*(1.0D0 - 6.0D0*beta*Xs*XA1)
       DN2 = 1.0D0/(1.0D0+Xs*Xs) + 2.0D0*Fb*(2.0D0 - 6.0D0*beta*Xs*XA1)
       DN3 = siper*(1.0D0+2.0D0*Fb) + XA1*DN2
       D02=D0/2.0D0
*      1.25992104989487030D0 = 2.0D0**0.33333333333333D0
       de1 = 1.33333333333333D0/(dens**2.33333333333333D0) 
       DGRADx = (Dx*Dxx + Dy*Dxy + Dz*Dxz)/DGRAD
       GRADXx = 2.0D0**0.33333333333333D0*
     >         (1/(dens*y)*DGRADx - de1*Dx*DGRAD)
       DGRADy = (Dx*Dxy + Dy*Dyy + Dz*Dyz)/DGRAD
       GRADXy = 2.0D0**0.33333333333333D0*
     >         (1/(dens*y)*DGRADy - de1*Dy*DGRAD)
       DGRADz = (Dx*Dxz + Dy*Dyz + Dz*Dzz)/DGRAD
       GRADXz = 2.0D0**0.33333333333333D0*
     >         (1/(dens*y)*DGRADz - de1*Dz*DGRAD)
       T1 = Dx/2.0D0*GRADXx
       T2 = Dy/2.0D0*GRADXy
       T3 = Dz/2.0D0*GRADXz
       DN4 = 6.0D0*beta*Fb*(T1+T2+T3)
       DN5 = 1.33333333333333D0*r43*r13*Xs*Xs
       TOT2 = DN5 - D02*DN1 + DN4*DN3
*       vxc = -2.0D0*beta*Fb/r43*TOT2
       vxc = -beta*Fb/r43*TOT2
       v = v0 + vxc 
*   check
*       return
*
c
c -----
c correlation Perdew
c local part ( equal to Vosko's)
c
  452   dens2=dens**2 
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
       if (IEXCH.eq.6) then
         v=v+vc
         return
        endif
c
c
       rs2=rs**2
       Cx1=0.002568D0+alf*rs+bet*rs2
       Cx2=1.D0+gam*rs+del*rs2+1.0D04*bet*rs**3
       C=0.001667D0+Cx1/Cx2
       Cx3=alf+2.D0*bet*rs
       Cx4=gam+2.D0*del*rs+3.0D4*bet*rs2
       dC=Cx3/Cx2-Cx1/Cx2**2*Cx4
       dC=-0.333333333333333D0*dC*const2/(y*dens)
c
c      phi=1.745*0.11*0.004235/C*DGRAD/dens**(7.D0/6.D0)
       phi=0.0008129082D0/C*DGRAD/dens**(1.1666666666666666D0)
       expo=exp(-phi)
       ex0=expo*C
       ec=ec+ex0*DGRAD2/(y*dens2)
*       e=e+ec
c
       D1=(2.D0-phi)*D0/dens
       phi2=phi**2
       D2=1.33333333333333333D0-3.666666666666666666D0*phi+
     >      1.166666666666666D0*phi2
       D2=D2*DGRAD2/dens2
       D3=phi*(phi-3.D0)*u0/(dens*DGRAD)
       D4=expo*DGRAD2/(y*dens)*(phi2-phi-1.D0)*dC
       vc=vc-1.0D0*(ex0/y*(D1-D2+D3)-D4)
       v=v+vc
c
       return

      end
c---------------------------------------------------------------------
c
      FUNCTION ASINH(X)
      implicit real*8 (a-h,o-z)
c
      t1=sqrt(X**2+1.D0)
      if (X.gt.0.0D0) then
       asinh=log(X+t1)
      else
       asinh=-log(-X+t1)
      endif
c
      return
      end
c-------------------------------------------------------------------

