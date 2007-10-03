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

         lap = dxx + dyy + dzz
 
          
           
            call closedpbe(dens,agrad,delgrad,lap,
     1           expbe,vxpbe,ecpbe,vcpbe)

  

            ex = expbe
            ec = ecpbe 
            v = vxpbe + vcpbe
 
C============DEBUG================
C
c           ex = 0.0D0
c           ec = 0.0D0
c           vxpbe = 0.0D0
c           vcpbe = 0.0D0
c           v = vxpbe + vcpbe    
C          
C==================================
 
           end

c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      subroutine closedpbe(rho,agrad,delgrad,lap,
     1           expbe,vxpbe,ecpbe,vcpbe)

      implicit real*8(a-h,o-z)
      parameter(thrd=1.d0/3.d0,thrd2=2.d0*thrd)
      parameter(pi32=29.608813203268075856503472999628d0)
      parameter(pi=3.1415926535897932384626433832795d0)
      parameter(alpha=1.91915829267751300662482032624669d0)
      parameter(thrd4=4.d0/3.d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
      parameter(thrdm=-thrd)
      parameter(sixthm=thrdm/2.d0)
      parameter(gama=0.03109069086965489503494086371273d0)
      parameter(bet=0.06672455060314922d0,delt=bet/gama)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c        IF(rho.lt.1.0D-18) THEN
c        expbe = 0.0D0
c        vxpbe = 0.0D0
c        ecpbe = 0.0D0
c        vcpbe = 0.0D0
c        ELSE
c        ENDIF
        
         rho2 = rho*rho
         rho13 = rho**thrd
         fk1 = pi32**thrd
         fk = fk1*rho13

C S = |grad(rho)|/(2*fk*rho)

         s1 = 2.0D0*fk*rho
          s = agrad / s1
           
C LDA exchange contribution:
C ex*rho = -(3/4Pi)*(e^2)*(3pi)^2/3*rho^1/3*rho
C ex*rho = -0.75*(3/Pi)^1/3*rho^4/3
C ex*rho = -ax*rho^4/3
 
         rho43 = rho**thrd4
         exlsd = ax*rho43        
                  
C In order to calculate the PBE contribution
C to exchange energy, we have to calculate the
C enhancement function Fx:
C Fx = 1+uk -(uk/(1+(um*s^2)/uk)

          s2 = s*s

C um/uk = ul
C P0 = 1 + (m*s^2)/uk

          p0 = 1.0d0 +ul*s2 
          fx = 1.0d0 + uk - uk/p0

C exchange pbe energy
C===================================================

          expbe = exlsd*fx

C====================================================

C Now the potential:
C dEx/drho = ax*rho^1/3*(4/3*Fx-t*(s^-1)(dF/ds)-(u-4/3*s^3)*(d/ds)*(s^-1*dF/ds))
C dEx/drho = ax*rho^1/3*(4/3*Fx-t*Fs-(u-4/3*s^3)*Fss)
C where u and t are:
C t = (2*kf)^-2*rho^-1*lap(rho)
C u = (2*kf)^-3(rho)^-2*grad(rho)*grad(|grad(rho|)

           s3 = s2*s
           twofk = 2.0D0*fk
           twofk2 = twofk*twofk
           twofk3 = twofk*twofk2
           v1 = twofk2*rho
           v = lap/t1
            
           u1 = twofk3*rho2   
           u = delgrad/u1

C Calculation of first and second derivatives
C Fs = (1/s)*(dF/ds)= (2*um)/p0^2
           
            P2 = P0*P0
            Fs = 2.0D0*um/P2

C Fss = dFs/ds
C Fss = -4.0D0*ul*s*Fs/P0

            F1 = -4.0D0*ul*s*Fs
            Fss = F1/P0


C Now we calculate the potential Vx
            
            vx1 = ax*rho13
            vx2 = thrd4*Fx
            vx3 = v*Fs  
            vx4 = (u-thrd4*s3)*Fss

C======================================================
                      
            vxpbe = vx1*(vx2-vx4-vx3)

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

                   CALL GCORc(rtrs,ec,ecrs)
C===================================================================
   
           eclsd = ec
           vclsd = eclsd - rs*thrd*ecrs
          
C===================================================================

C Now we have to calculate the H function in order to evaluate
C the GGA contribution to the correlation energy
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      H = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)

           PON = -eclsd / GAMA
           B = DELT / (DEXP(PON) - 1.0D0)
           B2 = B*B
           T4 = T2*T2
           RS2 = RS*RS
           RS3 = RS2*RS

            Q4 = 1.0D0+B*T2
            Q5 = 1.0D0+B*T2+B2*T4

             H = (BET/DELT)*DLOG(1.0D0+DELT*Q4*T2/Q5)

C So the correlation energy for pbe is:

C===============================================================

          ecpbe = eclsd + h

C==============================================================

C Now we have to calculate the potential contribution of GGA
      
      sk = DSQRT(4.0D0*fk/pi)
      twoks = 2.0D0*sk
      twoks2 = twoks*twoks
      twoks3 = twoks2*twoks
      UU = delgrad / (rho2*twoks3)
      VV = lap/(rho2*twoks2)
    
      T6 = T4*T2
      RSTHRD = RS/3.D0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/BET
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      hB = -BET*B*T6*(2.D0+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      hBT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      hRST = RSTHRD*T2*hBT*BEC*ECRS
      hT = 2.d0*BET*Q9/Q8
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      hTT = 4.D0*BET*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      COMM = COMM-UU*HTT-VV*HT

      DVC = COMM

C Then, the potential for PBE is:

C================================================================

      vcpbe = vclsd + dvc

C================================================================
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

