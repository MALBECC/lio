c------------------------------------------------------------------
c subroutine for evaluating exchange correlation density and
c potential , for non local density functionals
c Iexch, index should be greater than 4
c 4 : Exchange given by Perdew : Phys. Rev B 33 8800 (1986)
c     Correlation :     Perdew : Phys. Rev B 33 8822 (1986)
c
c 5 : Exchange given by Becke  : Phys. Rev A 38 3098 (1988)
c     Correlation :     Perdew : Phys. Rev B 33 8822 (1986) 
c to be used in concert with dnsgop subroutine, the one that gives
c gradient and second derivatives of the electronic density
c 19-1-93
c------------------------------------------------------------------
c
c     
      subroutine potgop(Iexch,densa,densb,adx,bdx,ady,bdy,adz,bdz,adxx,
     >           bdxx,adyy,bdyy,adzz,bdzz,adxy,bdxy,adyz,bdyz,adxz,bdxz,
     >           ex,ec,va,vb)
      implicit real*8 (a-h,o-z)
      external asinh
c
c data X alpha
      data const /-0.738558766382022447D0/
      data const3 /-0.930525736349100185D0/
c data Gunnarson-Lundvquist
      data const2 /0.620350490899400087D0 /
c  data Vosko et al
      data A1,b1,c1,x0,Q,A16 /0.03109205D0,3.72744D0,12.9352D0,
     >  -0.10498D0,6.15199066246304849D0, 0.0051820083D0 /
      data A2,b2,c2,x02,Q2,A26 /0.015546025D0,7.06042D0,18.0578D0,
     >  -0.32500D0,4.7309269D0,0.0025910042D0 /
      data A3,b3,c3,x03,Q3,A36 /.01554535D0,7.06042D0,18.0578D0,
     >     -0.32500,4.730926909D0,0.D0/
      data A4,b4,c4,x04,Q4,A46 /-0.016887D0,1.131071D0,13.0045D0,
     >     -0.0047584D0,7.123108759D0,0.D0/
*  data Becke
      data beta /0.0042D0/
c
c non local Perdew functional
c 
      data alf,bet,gam,del /0.023266D0,7.389D-06,8.723D0,0.472D0/
*
        dens = densa + densb
        dens2 = dens * dens
        dx = adx + bdx
        dy = ady + bdy
        dz = adz + bdz
        dxx = adxx + bdxx
        dyy = adyy + bdyy
        dzz = adzz + bdzz
        dxy = adxy + bdxy
        dyz = adyz + bdyz
        dxz = adxz + bdxz
c
        if (dens.eq.0.D0) then
          ex=0.0D0
          ec = 0.D0
          va=0.0D0
          vb=0.0D0
          return
        endif
c
        y=dens**0.333333333333333333D0
c
c
c-- Non local --------------------------- 
c
*
       DGRAD2 = Dx*Dx + Dy*Dy + Dz*Dz
       DGRAD = dsqrt(DGRAD2)
       u0 = (Dx*Dx*Dxx + 2.D0*Dx*Dy*Dxy + 2.D0*Dy*Dz*Dyz + 
     >       2.D0*Dx*Dz*Dxz + Dy*Dy*Dyy + Dz**2.0D0*Dzz)/DGRAD
       D0 = Dxx + Dyy + Dzz
c
       if(iexch.eq.5) goto 451
*
*   Perdew Exchange (Phys. Rev B, 33 8800)
*
*       fx = 1.D0/15.D0
       fx = 0.0666666666666667D0 
*
*  spin a calculations
       if (densa.eq.0.D0) then
          ea = 0.D0
          va = 0.D0
          ya = 0.D0
          goto 453
       endif
       ya = (2.D0*densa)**0.333333333333333333D0
       ckFa = 3.0936677D0*ya
       aDGRAD2 = 4.D0*(aDx*aDx + aDy*aDy + aDz*aDz)
       aDGRAD = dsqrt(aDGRAD2)
       u0a =8.D0*(aDx*aDx*aDxx + 2.D0*aDx*aDy*aDxy + 2.D0*aDy*aDz*aDyz +
     >      2.D0*aDx*aDz*aDxz + aDy*aDy*aDyy+ aDz*aDz*aDzz) / aDGRAD
       D0a = 2.D0*(aDxx + aDyy + aDzz)
       densa2 = 4.D0*densa*densa
       sa = aDGRAD/(4.D0*ckFa*densa)
c
       sa2 = sa*sa
       sa3 = sa*sa2
       sa4 = sa*sa3
       g0a = 1.D0 + 1.296D0*sa2 + 14.D0*sa4 + 0.2D0*sa3*sa3
       Fa = g0a**fx
       ea = densa/dens*const*Fa
*
       ta = D0a/(densa*8.D0*ckFa*ckFa)
       ua = u0a/((2.D0*ckFa)**3*densa2)
       g2a = 2.592D0*sa + 56.D0*sa3 + 1.2D0*sa2*sa3
       g3a = 2.592D0 + 56.D0*sa2 + 1.2D0*sa4
       g4a = 112.D0*sa + 4.8D0*sa3
       dFa = fx*Fa/g0a*g2a
       dsFa = fx*Fa/g0a*(-14.D0*fx*g3a*g2a/g0a + g4a)
       va = const*(1.33333333333D0*Fa - ta/sa*dFa -
     >      (ua - 1.3333333333D0*sa3)*dsFa)*ya
*
 
*  spin b calculations
 453   if (densb.eq.0.D0) then
          eb = 0.D0
          vb = 0.D0
          yb = 0.D0
          goto 454
       endif
          
       yb = (2.D0*densb)**0.333333333333333333D0
       ckFb = 3.0936677D0*yb
       bDGRAD2 = 4.D0*(bDx*bDx + bDy*bDy + bDz*bDz)
       bDGRAD = dsqrt(bDGRAD2)
       u0b =8.D0*(bDx*bDx*bDxx + 2.D0*bDx*bDy*bDxy + 2.D0*bDy*bDz*bDyz +
     >      2.D0*bDx*bDz*bDxz + bDy*bDy*bDyy+ bDz*bDz*bDzz) / bDGRAD
       D0b = 2.D0*(bDxx + bDyy + bDzz)
       densb2 = 4.D0*densb*densb
       sb = bDGRAD/(4.D0*ckFb*densb)
c
       sb2 = sb*sb
       sb3 = sb2*sb
       sb4 = sb3*sb
       g0b = 1.D0 + 1.296D0*sb2 + 14.D0*sb4 + 0.2D0*sb3*sb3
       Fb = g0b**fx
       eb = densb/dens*const*Fb
****
c
c
       tb = D0b/(densb*8.D0*ckFb*ckFb)
       ub = u0b/((2.D0*ckFb)**3*densb2)
c
       g2b = 2.592D0*sb + 56.D0*sb3 + 1.2D0*sb2*sb3
       g3b = 2.592D0 + 56.D0*sb2 + 1.2D0*sb4
       g4b = 112.D0*sb + 4.8D0*sb3
       dFb = fx*Fb/g0b*g2b
       dsFb = fx*Fb/g0b*(-14.D0*fx*g3b*g2b/g0b + g4b)
       vb = const*(1.33333333333D0*Fb - tb/sb*dFb -
     >      (ub - 1.3333333333D0*sb3)*dsFb)*yb
*
 454   ex = ya*ea + yb*eb
*
************   debug
*       return
*       e = 0.D0
*       va = 0.D0
*       vb = 0.D0
************
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
*          a spin   -   energy
 451   e = 0.D0
       ea = 0.D0
       eb = 0.D0
       va = 0.D0
       vb = 0.D0
       if ( densa.eq.0.D0) then
           ea = 0.0D0
           va = 0.0D0
           goto 455
       endif
       ya = densa**0.333333333333333333D0
       y43a = ya*densa
       e0a = const3*y43a
       aDGRAD2 = aDx*aDx + aDy*aDy + aDz*aDz
       aDGRAD = dsqrt(aDGRAD2) 
       Xsa = aDGRAD/y43a
       sipera = asinh(Xsa)
       DNa = 1.0D0 + 6.0D0*beta*Xsa*sipera
       eca = - beta*y43a*Xsa*Xsa/DNa
       ea = e0a + eca 
*
       v0a = 1.33333333333333D0*const3*ya
       Fba = 1.0D0/DNa
       XA1a = Xsa/dsqrt(1+Xsa*Xsa)
       DN1a = 1.0D0 + Fba*(1.0D0 - 6.0D0*beta*Xsa*XA1a)
       DN2a = 1.0D0/(1.0D0+Xsa*Xsa) +
     >        2.0D0*Fba*(2.0D0 - 6.0D0*beta*Xsa*XA1a)
       DN3a = sipera*(1.0D0+2.0D0*Fba) + XA1a*DN2a
       D0a = aDxx + aDyy + aDzz
       de1a = 1.33333333333333D0/(densa**2.33333333333333D0)
       aDGRADx = (aDx*aDxx + aDy*aDxy + aDz*aDxz)/aDGRAD
       aGRADXx = (1/(densa*ya)*aDGRADx - de1a*aDx*aDGRAD)
       aDGRADy = (aDx*aDxy + aDy*aDyy + aDz*aDyz)/aDGRAD
       aGRADXy = (1/(densa*ya)*aDGRADy - de1a*aDy*aDGRAD)
       aDGRADz = (aDx*aDxz + aDy*aDyz + aDz*aDzz)/aDGRAD
       aGRADXz = (1/(densa*ya)*aDGRADz - de1a*aDz*aDGRAD)
       T1a = aDx*aGRADXx
       T2a = aDy*aGRADXy
       T3a = aDz*aGRADXz
       DN4a = 6.0D0*beta*Fba*(T1a+T2a+T3a)
       DN5a = 1.33333333333333D0*y43a*ya*Xsa*Xsa
       TOT2a = DN5a - D0a*DN1a + DN4a*DN3a
       vxca = -beta*Fba/y43a*TOT2a
       va = v0a + vxca

*
*
*        b spin   -   energy
 455   if ( densb.eq.0.D0) then
           eb = 0.0D0
           vb = 0.0D0
           goto 456
       endif
       yb = densb**0.333333333333333333D0
       y43b = yb*densb
       e0b = const3*y43b
       bDGRAD2 = bDx*bDx + bDy*bDy + bDz*bDz
       bDGRAD = dsqrt(bDGRAD2)
       Xsb = bDGRAD/y43b
       siperb = asinh(Xsb)
       DNb = 1.0D0 + 6.0D0*beta*Xsb*siperb
       ecb = - beta*y43b*Xsb*Xsb/DNb
       eb = e0b + ecb
*
*        b spin   -   potential
       v0b = 1.33333333333333D0*const3*yb
       Fbb = 1.0D0/DNb
       XA1b = Xsb/dsqrt(1+Xsb*Xsb)
       DN1b = 1.0D0 + Fbb*(1.0D0 - 6.0D0*beta*Xsb*XA1b)
       DN2b = 1.0D0/(1.0D0+Xsb*Xsb) +
     >        2.0D0*Fbb*(2.0D0 - 6.0D0*beta*Xsb*XA1b)
       DN3b = siperb*(1.0D0+2.0D0*Fbb) + XA1b*DN2b
       D0b = bDxx + bDyy + bDzz
       de1b = 1.33333333333333D0/(densb**2.33333333333333D0)
       bDGRADx = (bDx*bDxx + bDy*bDxy + bDz*bDxz)/bDGRAD
       bGRADXx = (1/(densb*yb)*bDGRADx - de1b*bDx*bDGRAD)
       bDGRADy = (bDx*bDxy + bDy*bDyy + bDz*bDyz)/bDGRAD
       bGRADXy = (1/(densb*yb)*bDGRADy - de1b*bDy*bDGRAD)
       bDGRADz = (bDx*bDxz + bDy*bDyz + bDz*bDzz)/bDGRAD
       bGRADXz = (1/(densb*yb)*bDGRADz - de1b*bDz*bDGRAD)
       T1b = bDx*bGRADXx
       T2b = bDy*bGRADXy
       T3b = bDz*bGRADXz
       DN4b = 6.0D0*beta*Fbb*(T1b+T2b+T3b)
       DN5b = 1.33333333333333D0*y43b*yb*Xsb*Xsb
       TOT2b = DN5b - D0b*DN1b + DN4b*DN3b
       vxcb = -beta*Fbb/y43b*TOT2b
       vb = v0b + vxcb

*
  456  ex = (ea + eb)/dens
*
*   check
**
* 452      return
**
  452    continue 
*
c
c -----
c Perdew Correlation 
*
c Local part ( equal to Vosko's)
c
        rs=const2/y
        x1=dsqrt(rs)
        zeta = (densa - densb)/dens
        zeta4 = zeta**4
        gi = 1.125D0*( (1+zeta)**1.333333333333333D0 +
     >                 (1-zeta)**1.333333333333333D0 - 2)
        Xx=rs+b1*x1+c1
         Xx3=rs+b3*x1+c3
         Xx4=rs+b4*x1+c4
        Xxo=x0**2+b1*x0+c1
         Xxo3=x03**2+b3*x03+c3
         Xxo4=x04**2+b4*x04+c4
        t1=2.D0*x1+b1
         t13=2.D0*x1+b3
         t14=2.D0*x1+b4
        t2=log(Xx)
         t23=log(Xx3)
         t24=log(Xx4)
        t3=atan(Q/t1)
         t33=atan(Q3/t13)
         t34=atan(Q4/t14)
        t4=b1*x0/Xxo
         t43=b3*x03/Xxo3
         t44=b4*x04/Xxo4
*
       ecP=A1*(2.D0*log(x1)- t2+ 2.D0*b1/Q*t3 -t4*(2.D0*log(x1-x0)-
     >     t2+ 2.D0*(b1+2.D0*x0)/Q*t3))
       ecF=A3*(2.D0*log(x1)- t23+ 2.D0*b3/Q3*t33 -t43*(2.D0*log(x1-x03)-
     >     t23+ 2.D0*(b3+2.D0*x03)/Q3*t33))
       ecA=A4*(2.D0*log(x1)- t24+ 2.D0*b4/Q4*t34 -t44*(2.D0*log(x1-x04)-
     >     t24+ 2.D0*(b4+2.D0*x04)/Q4*t34))
        acca = 1.709920934161367D0*(ecF-ecP)/ecA -1
        ec = ecP + ecA*gi*(1 + acca*zeta4)
*
        t5 = t1/Xx
         t53 = t13/Xx3
         t54 = t14/Xx4
        t6 = 1/(t1**2 + Q**2)
         t63 = 1/(t13**2 + Q3**2)
         t64 = 1/(t14**2 + Q4**2)
       decP = A1*(2/x1 - t5 - 4.D0*b1*t6 - b1*x0/Xxo * (2/(x1-x0) -
     >         t5 - 4.D0*(+2.D0*x0+b1)*t6 ))
*       decP = A1*(2/x1 - t5 - 4.D0*b1*t6 - b1*x0/Xxo * (2/(x1-x0) -
*     >          t5 - 4.D0*(-2.D0*x0+b1)*t6 ))
       decF = A3*(2/x1 - t53 - 4.D0*b3*t63 - b3*x03/Xxo3 * (2/(x1-x03) -
     >         t53 - 4.D0*(2.D0*x03+b3)*t63 ))
       decA = A4*(2/x1 - t54 - 4.D0*b4*t64 - b4*x04/Xxo4 * (2/(x1-x04) -
     >         t54 - 4.D0*(2.D0*x04+b4)*t64 ))
*
        dacca = 1.709920934161367D0/ecA * ( decF - decP - (ecF-ecP)/
     >          ecA*decA )
        dgi = 1.5D0*((1+zeta)**0.3333333333333333D0 -
     >              (1-zeta)**0.3333333333333333D0 )
        dazeta = 1/dens*(1-zeta)
        dbzeta = -1/dens*(1+zeta)
        term1 = -x1/6.D0*(decP + decA*gi*(1+acca*zeta4) +
     >        ecA*gi*dacca*zeta4)
        term2 = ecA*(dgi*(1+acca*zeta4) +
     >        4*gi*acca*zeta**3)
        term3 = term1 + term2*dazeta*dens
        term4 = term1 + term2*dbzeta*dens
        daec = ec + term3
        dbec = ec + term4
*
*        e = e + ec
        va = va + daec
        vb = vb + dbec
*
*  Non-local part 
*
      zeta = (densa - densb)/dens
      di = 2**0.333333333333333D0*((0.5D0*(1+zeta))**1.666666666666666D0
     >     + (0.5D0*(1-zeta))**1.666666666666666D0 ) **0.5D0
       rs2 = rs**2
       Cx1 = 0.002568D0 + alf*rs + bet*rs2
       Cx2 = 1.D0 + gam*rs + del*rs2 + 1.0D04*bet*rs**3
       C = 0.001667D0 + Cx1/Cx2
       Cx3 = alf + 2.D0*bet*rs
       Cx4 = gam + 2.D0*del*rs + 3.0D4*bet*rs2
       dC = Cx3/Cx2 - Cx1/Cx2**2*Cx4
       dC = -0.333333333333333D0*dC*const2/(y*dens)
c
c      phi=1.745*0.11*0.004235/C*DGRAD/dens**(7.D0/6.D0)
       phi = 0.0008129082D0/C*DGRAD/dens**(1.1666666666666666D0)
       expo = exp(-phi)
       ex0 = expo*C
       ecnl = ex0*DGRAD2/(y*dens2)/di
*       e = e + ecnl
       ec = ec + ecnl 
c
       D1 = (2.D0 - phi)*D0/dens
       phi2 = phi**2
       D2 = 1.33333333333333333D0 - 3.666666666666666666D0*phi +
     >      1.166666666666666D0*phi2
       D2 = D2*DGRAD2/dens2
       D3 = phi*(phi - 3.D0)*u0/(dens*DGRAD)
       D4 = expo*DGRAD2/(y*dens)*(phi2 - phi - 1.D0)*dC
      D5 = 0.83333333333333333D0/(di*di*dens**3.666666666666666666D0)
     >  *(densa**0.666666666666666666D0 - densb**0.666666666666666666D0)
      D6a = 1.587401052D0*(1.D0 - phi)*densb*DGRAD2
      D6b = 1.587401052D0*(1.D0 - phi)*densa*DGRAD2
      D7 = 1.587401052D0*(2.D0 - phi)*dens
      D8a = bDx*Dx + bDy*Dy + bDz*Dz
      D8b = aDx*Dx + aDy*Dy + aDz*Dz
*
       vca =  - (ex0/y*(D1 - D2 + D3 - D5*(D6a + D7*D8a)) - D4) / di
       vcb =  - (ex0/y*(D1 - D2 + D3 + D5*(D6b + D7*D8b)) - D4) / di
*
       va = va + vca
       vb = vb + vcb
c
       return

      end
