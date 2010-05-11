
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
c 7 : BLYP, Exchange: given by Becke 
c           Correlation: given by LYP: PRB 37 785 (1988)
c 8: Exchange given by Perdew : Phys. Rev B 33 8800 (1986)
c       Correlation: given by LYP: PRB 37 785 (1988)
c to be used in concert with dnsg subroutine, the one that gives
c gradient and second derivatives of the electronic density
c 19-1-93, last version: 25/11/97.
c 9: PBE
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
c  data for LYP functional (correlation)
      data alyp,blyp,clyp,dlyp /0.04918d0,0.132d0,0.2533d0,0.349d0/
      data dlyp3,clyp3,cf /0.116333333D0,0.0844333333D0,2.871234D0/
c
c definition of pi:
       pi=acos(-1.d0)
c
      IF(IEXCH.GE.4.AND.IEXCH.LE.9) THEN
      CONTINUE
      ELSE
      WRITE(*,*) 'OPCION DE FUNCIONAL NO VALIDA'
      STOP
      ENDIF

c        if (dens.eq.0) then
      if (dens.lt.1D-13) then
         v=0.0D0
         ex=0.0D0
         ec = 0.D0
         return
        endif
c
        y=dens**0.333333333333333333D0
        v=0
c
c
c-- Non local --------------------------- 
c
C  Modified By D. Bikiel-------------------

        IF(IEXCH.EQ.9) GOTO 460
c-------------------------------------------
c  Exchange Perdew, Correlation LYP
c
       if (iexch.eq.8) then
c Perdew (Phys. Rev B, 33 8800)
c results checked for energies against paper
       ckF=3.0936677D0*y
       DGRAD2=Dx**2+Dy**2+Dz**2
       DGRAD=sqrt(DGRAD2)
       u0=(Dx**2*Dxx+ 2.D0*Dx*Dy*Dxy+ 2.D0*Dy*Dz*Dyz+ 2.D0*Dx*Dz*Dxz+
     >   Dy**2.0D0*Dyy+ Dz**2.0D0*Dzz)/DGRAD
       D0=Dxx+Dyy+Dzz
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
c 
c Correlation BLYP
         rom13=dens**(-0.3333333333D0)
          rom53=dens**(1.666666666666D0)
          ecro=dexp(-clyp*rom13)
          f1=1.d0/(1.d0+dlyp*rom13)
c         cf=3.d0/10.d0*(3.d0*pi**2)**(2.d0/3.d0)
          tw=1.d0/8.d0*(DGRAD2/dens-D0)
          term=(tw/9.d0+D0/18.d0)-2.d0*tw+cf*(rom53)
          term=dens + blyp*(rom13**2)*ecro*term
c energy:
          ec=-alyp*f1*term/dens
c
          h1=ecro/rom53
          g1=f1*h1
          tm1=dlyp3*(rom13/dens)
          fp1=tm1*(f1**2)
          tm2=-1.666666666D0+clyp3*rom13
          hp1=h1*tm2/dens
          gp1=fp1*h1+hp1*f1
          fp2=tm1*2.0D0*f1*(fp1-0.6666666666D0*f1/dens)
          tm3=1.6666666666D0-clyp3*1.3333333333D0*rom13
          hp2=hp1*tm2/dens+h1*tm3/dens**2
          gp2=fp2*h1+2.0D0*fp1*hp1+hp2*f1

c
          term3=-alyp*(fp1*dens+f1)-alyp*blyp*cf*(gp1*dens+
     >           8.d0/3.d0*g1)*rom53
          term4=(gp2*dens*DGRAD2+gp1*(3.d0*DGRAD2+2.d0*dens*D0)+
     >           4.d0*g1*D0)*alyp*blyp/4.d0
          term5=(3.d0*gp2*dens*DGRAD2+gp1*(5.d0*DGRAD2+6.d0*dens*D0)+
     >           4.d0*g1*D0)*alyp*blyp/72.d0
c potential:
          vc=term3-term4-term5
          v = v + vc
          return
        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
c       write(*,*) t,u
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
c       write(*,*) Dxx,Dyy,Dzz
*   check
*       return
*
c------------------------------------------
c correlation LYP (option: 7)
c
       if (iexch.eq.7) then
          rom13=dens**(-0.3333333333D0)
          rom53=dens**(1.666666666666D0)
          ecro=dexp(-clyp*rom13)
          f1=1.d0/(1.d0+dlyp*rom13)
c         cf=3.d0/10.d0*(3.d0*pi**2)**(2.d0/3.d0)
          tw=1.d0/8.d0*(DGRAD2/dens-D0)
          term=(tw/9.d0+D0/18.d0)-2.d0*tw+cf*(rom53)
          term=dens + blyp*(rom13**2)*ecro*term
c energy:
          ec=-alyp*f1*term/dens
c
          h1=ecro/rom53
          g1=f1*h1
          tm1=dlyp3*(rom13/dens)
          fp1=tm1*(f1**2)
          tm2=-1.666666666D0+clyp3*rom13
          hp1=h1*tm2/dens
          gp1=fp1*h1+hp1*f1
          fp2=tm1*2.0D0*f1*(fp1-0.6666666666D0*f1/dens)
          tm3=1.6666666666D0-clyp3*1.3333333333D0*rom13
          hp2=hp1*tm2/dens+h1*tm3/dens**2
          gp2=fp2*h1+2.0D0*fp1*hp1+hp2*f1         

c
          term3=-alyp*(fp1*dens+f1)-alyp*blyp*cf*(gp1*dens+
     >           8.d0/3.d0*g1)*rom53
          term4=(gp2*dens*DGRAD2+gp1*(3.d0*DGRAD2+2.d0*dens*D0)+
     >           4.d0*g1*D0)*alyp*blyp/4.d0
          term5=(3.d0*gp2*dens*DGRAD2+gp1*(5.d0*DGRAD2+6.d0*dens*D0)+
     >           4.d0*g1*D0)*alyp*blyp/72.d0
c potential:
          vc=term3-term4-term5         
          v = v + vc
          return
        endif
c -----------------------------------------
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
c        write(*,*) 2.D0*log(x1),t2,t3,t4
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
c       write(*,*) ec,v,D1,D2,D3
c
       return

C==============================================================================
C PBE and PW91 functional incorporated from Burke's code. Minimal modification
C to run under this DFT program by D. Bikiel.
C==============================================================================
C##############################################################################

C  Density magnitudes

460      grad2 = dx*dx + dy*dy + dz*dz
         agrad = DSQRT(grad2)

         dgrad1 = dx*dx*dxx / agrad
         dgrad2 = dy*dy*dyy / agrad
         dgrad3 = dz*dz*dzz / agrad
         dgrad4 = 2.0D0*dx*dy*dxy / agrad
         dgrad5 = 2.0D0*dx*dz*dxz / agrad
         dgrad6 = 2.0D0*dy*dz*dyz / agrad
          
         delgrad = dgrad1+dgrad2+dgrad3+dgrad4+dgrad5+dgrad6

         rlap = dxx + dyy + dzz
 
           
            call closedpbe(dens,agrad,delgrad,rlap,
     1           expbe,vxpbe,ecpbe,vcpbe)
           

  

            ex = expbe
            ec = ecpbe 
            v = vxpbe + vcpbe

c        write(*,*)'ex,vx,ec,vc',ex,vxpbe,ec,vcpbe

C============DEBUG================
C
c           ex = 0.0D0
c           ec = 0.0D0
c           vxpbe = 0.0D0
c           vcpbe = 0.0D0
C          
c           v = vxpbe + vcpbe    
C==================================
 
           end

c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      subroutine closedpbe(rho,agrad,delgrad,rlap,
     1           expbe,vxpbe,ecpbe,vcpbe)

      implicit real*8(a-h,o-z)
      parameter(thrd=1.d0/3.d0)
      parameter(pi32=29.608813203268075856503472999628d0)
      parameter(pi=3.1415926535897932384626433832795d0)
      parameter(thrd4=4.d0/3.d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
      parameter(gama=0.03109069086965489503494086371273d0)
      parameter(bet=0.06672455060314922d0,delt=bet/gama)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
       
      IF(rho.LT.2.0D-18) THEN
      expbe = 0.0D0
      vxpbe = 0.0D0
      ecpbe = 0.0D0
      vcpbe = 0.0D0
      RETURN
      ENDIF

 
         rho2 = rho*rho
         rho13 = rho**thrd

         fk1 = pi32**thrd
         fk = fk1*rho13

C S = |grad(rho)|/(2*fk*rho)

          twofk = 2.0D0*fk
          twofk2 = twofk*twofk
          twofk3 = twofk2*twofk

          s = agrad / (twofk*rho)
          s2 = s*s  
          s3 = s2*s

C LDA exchange contribution:
C ex*rho ==> energy, we will calculate ex ==> energy density
C ex*rho = -(3/4Pi)*(e^2)*(3pi)^2/3*rho^1/3*rho
C ex*rho = -0.75*(3/Pi)^1/3*rho^4/3
C ex*rho = ax*rho^4/3
 
         exlda = ax*rho13        
                  
C In order to calculate the PBE contribution
C to exchange energy, we have to calculate the
C enhancement function Fx:
C Fx = 1+uk -(uk/(1+(um*s^2)/uk)
C um/uk = ul
C P0 = 1 + (um*s^2)/uk

          p0 = 1.0d0 + ul*s2 
          fxpbe = 1.0d0 + uk - uk/p0

C exchange pbe energy
C===================================================

          expbe = exlda*fxpbe
 
c        write(*,*)'exlda,fxpbe,expbe',exlda,fxpbe,expbe

C====================================================

C Now the potential:
C dEx/drho = ax*rho^1/3*(4/3*Fx-t*(s^-1)(dF/ds)-(u-4/3*s^3)*(d/ds)*(s^-1*dF/ds))
C dEx/drho = ax*rho^1/3*(4/3*Fx-t*Fs-(u-4/3*s^3)*Fss)
C where u and t are:
C t = (2*kf)^-2*rho^-1*lap(rho)
C u = (2*kf)^-3(rho)^-2*grad(rho)*grad(|grad(rho|)
          
c          write(*,*) lap    
         
           v = rlap/(twofk2*rho)
           u = delgrad/(twofk3*rho2)   

C Calculation of first and second derivatives
C Fs = (1/s)*(dF/ds)= (2*um)/p0^2
           
            P2 = P0*P0
            Fs = 2.0D0*um/P2

C Fss = dFs/ds
C Fss = -4.0D0*ul*s*Fs/P0

            F1 = -4.0D0*ul*s*Fs
            Fss = F1/P0


C Now we calculate the potential Vx
            
            vx2 = thrd4*Fxpbe
            vx3 = v*Fs  
            vx4 = (u-thrd4*s3)*Fss

C======================================================
                      
            vxpbe = exlda*(vx2-vx4-vx3)
c           write(*,*)'v,u,fxpbe,Fs,Fss',v,u,fxpbe,Fs,Fss  
C======================================================

C---------------------------------------------------
C Now we need to calculate the Correlation contribution
C to the energy
C ecpbe = eclsd*rho + h*rho
C first we calculate the lsd contribution to the correlation energy
C we will use the subroutine GCOR.
C We need only the  rs (seitz radius) rs = (3/4pi*rho)^1/3
 
      pirho = 4.0D0*Pi*rho
      rs = (3.0D0/pirho)**thrd
      rtrs = DSQRT(rs)

      sk = DSQRT(4.0D0*fk/pi)
      twoks = 2.0D0*sk

      t = agrad / (twoks*rho)
      t2 = t*t

      twoks2 = twoks*twoks
      twoks3 = twoks2*twoks

      UU = delgrad / (rho2*twoks3)
      VV = rlap/(rho*twoks2)


                   CALL GCORc(rtrs,ec,eurs)
C===================================================================
   
           eclda = ec
           ecrs = eurs
           vclda = eclda - rs*thrd*ecrs
c          vclda = 2.0D0*(eclda - rs*thrd*ecrs)
          
C===================================================================

C Now we have to calculate the H function in order to evaluate
C the GGA contribution to the correlation energy

           PON = -1.0D0*ec / GAMA
           B = DELT / (DEXP(PON) - 1.0D0)
           B2 = B*B
           T4 = T2*T2
           RS2 = RS*RS
           RS3 = RS2*RS

            Q4 = 1.0D0+B*T2
            Q5 = 1.0D0+B*T2+B2*T4

             H = (BET/DELT)*DLOG(1.0D0+DELT*Q4*T2/Q5)
c             write(*,*) DEXP(PON),T2

C So the correlation energy for pbe is:

C===============================================================

          ecpbe = eclda + h
c         write(*,*)'ecpbe,eclda,h',ecpbe,eclda,h
C==============================================================

C Now we have to calculate the potential contribution of GGA
      
      T6 = T4*T2
      RSTHRD = RS/3.D0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/BET
      BEC = B2*FAC/BET
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      hB = -BET*B*T6*(2.D0+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      hBT = 2.D0*BET*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      hRST = RSTHRD*T2*hBT*BEC*ECRS
      hT = 2.d0*BET*Q9/Q8
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      hTT = 4.D0*BET*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
 
c     write(*,*)'H,HRS,HRST,HT,T2,HTT,HBT',H,HRS,HRST,HT,T2,HTT,HBT 
       
      COMM = COMM-UU*HTT-VV*HT

      DVC = COMM
c     DVC = 2.0D0*COMM

C Then, the potential for PBE is:

C================================================================

      vcpbe = vclda + dvc
C================================================================


C===========DEBUG=================================
         
c         write(*,*)'rho,agrad,delgrad,lap',rho,agrad,delgrad,lap
c        write(*,*)'ex,vx',expbe,vxpbe
c        write(*,*)'ec,vc',ecpbe,vcpbe
c         write(*,*) expbe,Q4,H,eclda,ecpbe,vxpbe
c=====================================================

      END

C======================================================================
C#####################################################################
C=======================================================================
      SUBROUTINE GCORc(rtrs,GG,GRRS)

c slimmed down version of GCOR used in PW91 routines, to interpolate
c LSD correlation energy, as given by (10) of
c J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
c K. Burke, May 11, 1996.
C
C This subroutine calculates ec and its first derivative

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(A=0.0310907D0)
      PARAMETER(A1=0.21370D0)
      PARAMETER(B1=7.5957D0)
      PARAMETER(B2=3.5876D0)
      PARAMETER(B3=1.6382D0)
      PARAMETER(B4=0.49294D0)

      Q0 = -2.D0*A*(1.D0+A1*rtrs*rtrs)
      Q1 = 2.D0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2.D0*B2+rtrs*(3.D0*B3+4.D0*B4*rtrs))
      GRRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.d0+Q1))
      RETURN
      END

C=======================================================================
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

