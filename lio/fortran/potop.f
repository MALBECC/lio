c------------------------------------------------------------------
c OPEN SHELL VERSION
c subroutine for evaluating exchange correlation density and
c potential , depends of Iexch, index that says which potential
c will be used
c 1 X alpha
c 2 Gunnarson-Lundqvist
c 3 Vosko, Wilk and Nusair
c Exchange part is in all cases the same, for the time being LD
c Self interaction corrections could be used in correlation part
c------------------------------------------------------------------
c
c     
      subroutine potop(Iexch,densa,densb,ex,ec,va,vb)
      implicit real*8 (a-h,o-z)
c
c data X alpha
      data const,constv /-0.930525736349100185D0,-1.24070098179880017D0/

c data Gunnarson-Lundvquist
      data const2 /0.620350490899400087D0 /
c  data Vosko et al
      data A1,b1,c1,x0,Q,A16 /0.03109205D0,3.72744D0,12.9352D0,
     >     -0.10498D0,6.15199066246304849D0, 0.0051820083D0 /
      data A2,b2,c2,x02,Q2,A26 /0.015546025D0,7.06042D0,18.0578D0,
     >     -0.32500D0,4.7309269D0,0.0025910042D0 /
      data A3,b3,c3,x03,Q3,A36 /.01554535D0,7.06042D0,18.0578D0,
     >     -0.32500,4.730926909D0,0.D0/
      data A4,b4,c4,x04,Q4,A46 /-0.016887D0,1.131071D0,13.0045D0,
     >     -0.0047584D0,7.123108759D0,0.D0/
c
      dens=densa+densb
c
        if (dens.eq.0.0) then
         ex=0.0D0
         ec= 0.D0
         va=0.0D0
         vb=0.0D0
         return
        endif
c       
        ya=densa**0.333333333333333333D0 
        yb=densb**0.333333333333333333D0 
        y=densa*ya + densb*yb
        e0=const*y/dens
c       
        v0a=constv*ya
        v0b=constv*yb
c
        goto (10,20,30) Iexch
*
  10    ex = e0
        ec = 0.D0
        va=v0a
        vb=v0b
        return
c
  20    continue 
*
  30    ex = e0 
        y2 = dens**0.333333333333333333D0
        rs=const2/y2
        x1=sqrt(rs)
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
*        e = e0 + ec
        va = v0a + daec
        vb = v0b + dbec 
*
        return
      end
