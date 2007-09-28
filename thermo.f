c This little program is intended to be used in the cases
c in which the user wants to get the thermodynamical
c properties of a molecule and did not run the normal modes
c analysis with that option.
c
c It allows also  to include symmetry considerations (i.e CO
c entropy at 0K)
c examples : H2O,CO2 2
c homoatomic diatomic molecules also 2
c
c complex molecules, generally 1
c
c It's also convenient in cases in which information at several
c temperatures is needed.
c
c It reads all the information it needs from the
c file.vib file (obtained from a normal modes run)
c
c Formulas in Hill, 'An introduction to Statistical Thermodynamics'
c Dover, New York
c
c D.Estrin, Bs As Enero 1994
c----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      INCLUDE 'param'
      dimension r(nt,3),Pm(nt),vib(3*nt),Ei(3),Xi(20)
      dimension xmass(216)
      character*15 name4,name5,name1
c
      dimension Iz(nt),isotop(54)
      common /masses/ xmass
c
c reads input file (name of the original file)
      write(*,*) ' Name of the input file ?'
      read(*,100) name1
c
      ikk=1
      do 1 while (name1(ikk:ikk).ne.' ')
       ikk=ikk+1
 1    continue
c
      ikk=ikk-1
c
          name4=name1(1:ikk)//'.thermo'
          name5=name1(1:ikk)//'.vib'
c
      write(*,*) 'TEMPERATURE  and SIGMA'
      read(*,*) TEMP,sigma
   
      open(unit=2,file=name5)
      open(unit=44,file=name4)
      read(2,*) 
      read(2,*) natom
c
      do i = 1,54
        isotop(i) = 1
      enddo
c
      do n=1,natom
       read(2,*) Iz(n),r(n,1),r(n,2),r(n,3)
       indmass = (Iz(n)-1)*4 + isotop(Iz(n))
       Pm(n) = xmass(indmass)
      enddo
c
      Xcm=0.0D0
      Ycm=0.0D0
      Zcm=0.0D0
      Rcm=0.0D0
c
      do i=1,natom
       write(*,*) i,Pm(i)
       Rcm=Rcm+Pm(i)
       Xcm=Xcm+Pm(i)*r(i,1)
       Ycm=Ycm+Pm(i)*r(i,2)
       Zcm=Zcm+Pm(i)*r(i,3)
      enddo
c
      Xcm=Xcm/Rcm
      Ycm=Ycm/Rcm
      Zcm=Zcm/Rcm
c
      do i=1,10
       Xi(i)=0.0D0
      enddo
c
      do i=1,natom
       x1=r(i,1)-Xcm
       x2=x1**2
       y1=r(i,2)-Ycm
       y2=y1**2
       z1=r(i,3)-Zcm
       z2=z1**2
       r2=x2+y2+z2
c
       Xi(1)=Xi(1)+pm(i)*(r2-x2)
       Xi(4)=Xi(4)+pm(i)*(r2-y2)
       Xi(6)=Xi(6)+pm(i)*(r2-z2)
c
       Xi(2)=Xi(2)-pm(i)*x1*y1
       Xi(3)=Xi(3)-pm(i)*x1*z1
       Xi(5)=Xi(5)-pm(i)*y1*z1
c
      enddo
c
c ESSL OPTION ----------------------------------------
c#ifdef essl
c     call DSPEV(1,Xi,Ei,Wi,3,3,Xi(7),6)
c#endif
c
c EISPACK OPTION -----------------------------------------------
c#ifdef pack
c
       call dspev('N','L',3,Xi,Ei,Wi,3,Xi(7),info)
c#endif
c-----------------------------------------------------------
c Ei(1), Ei(2) and Ei(3) three principal inertia moments
c if molecule is linear one should be zero and the other two equal
c in atoms 3 of them are 0
      if(abs(Ei(1)).lt.1.0D-06) then
        if (abs(Ei(2)).lt.1.0D-06) then
         ntrasrot =3
        else
         ntrasrot = 5
        endif
      else
        ntrasrot = 6
      endif
      nnull = nt3 - nat3 + ntrasrot +1
      itemp = ntrasrot + 1
c
*
 
*  Conversion Factors :
*                      const1 converts frequency to cm**-1
*                      const2 converts intensity to km/mol
*                      const3 converts z.p. energy to kcal/mol
      const1 = 5141.1182009496200901D0
      const2 = 150.89915D0
      const3 = 1.430D-3
c
      const4= 1.43*2.0D0*4.186/(8.31451*TEMP)
      const5=TEMP*0.011545543
      const6=-1.164856775
      const7=2.072364943
c
c ALL THIS INFORMATION CORRESPONDS TO THE FOLLOWING MODEL:
c 1) RIGID ROTOR
c 2) HARMONIC UNCOUPLED OSCILLATORS
c 3) NO INTERNAL ROTATION
c 4) NO SYMMETRY FACTORS
c
c ENTROPY IN  IN CAL/K/MOL - ENTHAPLY CAL/MOL - HEAT CAPACITY
c CAL/K/MOL
c
c
c Translational entropy and enthalpy
c enthalpy : energy + RT, Pressure considered 1 atm
      Prs=1.0D0
      St=const6+1.50D0*log(Rcm)+2.50D0*log(TEMP)-log(Prs)
      St=St*8.31451/4.186
c
      Et=2.50D0*8.31451/4.186*TEMP
      Ct=1.50D0*8.31451/4.186
c
c Rotational entropy and energy
c
c      write(*,*) 'Inertia moments',Ei(1),Ei(2),Ei(3)
       if (abs(Ei(1)).lt.1.0D-06) then
        if (abs(Ei(2)).gt.1.0D-04) then
c linear molecules
         Er=8.31451/4.186*TEMP
c        Sr=const7+0.50D0*log(const5**3*Ei(2)*Ei(3))
         Sr=log(const5*Ei(2)/sigma)+1.0D0
         Cr=8.31451/4.186
        else
c atoms
         Er=0.0D0
         Sr=0.0D0
         Cr=0.0D0
        endif
       else
c non-linear molecules
       Er=1.50D0*8.31451/4.186*TEMP
       Sr=const7+0.50D0*log(const5**3*Ei(1)*Ei(2)*Ei(3))-log(sigma)
       Cr=1.50D0*8.31451/4.186
       endif
c
      Sr=Sr*8.31451/4.186
c Vibrational entropy
c
      read(2,*) 
      read(2,*)
      read(2,*)
      read(2,*)
c
      zpecm = 0.0D0
      do i=itemp,3*natom
       read(2,*)
       read(2,*) ki,freqcm, xintensity      
       vib(i)=freqcm*const4
       zpecm = zpecm + freqcm
        do n=1,natom
         read(2,*) x,y,z
        enddo
       read(2,*) 
      enddo
c Zero Point energy (in kcal/mol)
      zpecm = zpecm * const3

      Sv=0.0D0
      Ev=0.0D0
      Cv=0.0D0
       do i=itemp,3*natom
      if (vib(i).gt.100.D0) then
       rterm=0.0D0
       rterm2=0.0D0
      else
       rt1=exp(vib(i))
       rterm=vib(i)/(rt1-1.0D0)
       rterm2=rt1/(rt1-1.0D0)**2
      endif
      Sv=Sv+rterm-log(1.0D0-exp(-vib(i)))
c first expression includes zero point energy
c     Ev=Ev+0.50D0*vib(i)+rterm
      Ev=Ev+rterm
      Cv=Cv+vib(i)**2*rterm2
      enddo
      Sv=Sv*8.31451/4.186
      Ev=Ev*8.31451/4.186*TEMP
      Cv=Cv*8.31451/4.186
c
      write(44,700) TEMP,Prs
      write(44,*)
      write(44,701)
      write(44,702)
c
      write(44,703) Et,Ct,St
      write(44,704) Er,Cr,Sr
      write(44,705) Ev,Cv,Sv
      write(44,706) Et+Er+Ev,Ct+Cr+Cv,St+Sr+Sv
c
      write(44,*)
      write(44,707) zpecm
c
 100  format (A8)
 700  format(/,'THERMODYNAMIC PROPERTIES AT TEMPERATURE = ',f6.1,'K',3x,
     > 'PRESSURE= ',f5.2,' ATM')
 701  format(/,'CONTRIBUTION',5x,'ENTHALPY',5x,'HEAT CAPACITY',
     >  5x,'ENTROPY')
 702  format(17x,'CAL/MOL',6x,'CAL/MOL/K',9x,'CAL/MOL/K')
 703  format('TRANSLATIONAL',3x,f10.4,4x,f8.4,10x,f8.4)
 704  format('ROTATIONAL',6x,f10.4,4x,f8.4,10x,f8.4)
 705  format('VIBRATIONAL',5x,f10.4,4x,f8.4,10x,f8.4)
 706  format('     TOTAL ',5x,f10.4,4x,f8.4,10x,f8.4)
 707  format('ZERO POINT ENERGY (KCAL/MOL) = ',f8.3)
c
      close(44)

      end
c-------------------------------------------------------------------------
      BLOCK DATA
      implicit real*8 (a-h,o-z)
      common /masses/ xmass(216)
*     Atomic masses (u.m.a.) of most common isotopes
      data xmass /
*      H-1             H-2             H-3
     >  1.007825037D0,  2.014101787D0,  3.016049286D0,  0.0,
*      He-4            He-3
     >  4.00260325D0,   3.016029297D0,  0.0,            0.0,
*      Li-7            Li-6
     >  7.0160045D0,    6.0151232D0,    0.0,            0.0,
*      Be-9
     >  9.0121825D0,    0.0,            0.0,            0.0,
*       B-11            B-10
     > 11.0093053D0,   10.0129380D0,    0.0,            0.0,
*       C-12            C-13
     > 12.000000000D0, 13.003354839D0,  0.0,            0.0,
*       N-14            N-15
     > 14.003074008D0, 15.000108978D0,  0.0,            0.0,
*       O-16            O-18            O-17
     > 15.99491464D0,  17.99915939D0,  16.9991306D0,    0.0,
*       F-19
     > 18.99840325D0,   0.0,            0.0,            0.0,
*      Ne-20           Ne-22           Ne-21
     > 19.9924391D0,   21.9913837D0,   20.9938453D0,    0.0,
*      Na-23
     > 22.9897697D0,    0.0,            0.0,            0.0,
*      Mg-24           Mg-26           Mg-25
     > 23.9850450D0,   25.9825954D0,   24.9858392D0,    0.0,
*      AL-27
     > 26.9815413D0,    0.0,            0.0,            0.0,
*      Si-28           Si-29           Si-30
     > 27.9769284D0,   28.9764964D0,   29.9737717D0,    0.0,
*       P-31
     > 30.9737634D0,    0.0,            0.0,            0.0,
*       S-32            S-34            S-33            S-36
     > 31.9720718D0,   33.96786774D0,  32.9714591D0,   35.9670790D0,
*      Cl-35           Cl-37
     > 34.968852729D0, 36.965902624D0,  0.0,            0.0,
*      Ar-40           Ar-36           Ar-38
     > 39.9623831D0,   35.967545605D0, 37.9627322D0,    0.0,
*       K-39            K-41            K-40
     > 38.9637079D0,   40.9618254D0,   39.9639988D0,    0.0,
*      Ca-40           Ca-44           Ca-42           Ca-48
     > 39.9625907D0,   43.9554848D0,   41.9586218D0,   47.952532D0,
*      Sc-45
     > 44.9559136D0,    0.0,            0.0,            0.0,
*      Ti-48           Ti-46           Ti-47           Ti-49
     > 47.9479467D0,   45.9526327D0,   46.9517649D0,   48.9478705D0,
*       V-51            V-50
     > 50.9439625D0,   49.9471613D0,    0.0,            0.0,
*      Cr-52           Cr-53           Cr-50           Cr-54
     > 51.9405097D0,   52.9406510D0,   49.9460463D0,   53.9388822D0,
*      Mn-55
     > 54.9380463D0,    0.0,            0.0,            0.0,
*      Fe-56           Fe-54           Fe-57           Fe-58
     > 55.9349393D0,   53.9396121D0,   56.9353957D0,   57.9332778D0,
*      Co-59
     > 58.9331978D0,    0.0,            0.0,            0.0,
*      Ni-58           Ni-60           Ni-62           Ni-61
     > 57.9353471D0,   59.9307890D0,   61.9283464D0,   60.9310586D0,
*      Cu-63           Cu-65
     > 62.9295992D0,   64.9277924D0,    0.0,            0.0,
*      Zn-64           Zn-66           Zn-68           Zn-67
     > 63.9291454D0,   65.9260352D0,   67.9248458D0,   66.9271289D0,
*      Ga-69           Ga-71
     > 68.9255809D0,   70.9247006D0,    0.0,            0.0,
*      Ge-74           Ge-72           Ge-70           Ge-73
     > 73.9211788D0,   71.9220800D0,   69.9242498D0,   72.9234639D0,
*      As-75
     > 74.9215955D0,    0.0,            0.0,            0.0,
*      Se-80           Se-78           Se-82           Se-76
     > 79.9165205D0,   77.9173040D0,   81.916709D0,    75.9192066D0,
*      Br-79           Br-81
     > 78.9183361D0,   80.916290D0,     0.0,            0.0,
*      Kr-84           Kr-86           Kr-82           Kr-83
     > 83.9115064D0,   85.910614D0,    81.913483D0,    82.914134D0,
*      Rb-85
     > 84.9117D0,      0.0,             0.0,             0.0,
*      Sr-88           Sr-84           Sr-86           Sr-87
     > 87.9056D0,      83.9134d0,      85.9094d0,      86.9089d0,
*      Y-89
     > 88.9054D0,      0.0,             0.0,             0.0,
*      Zr-90           Zr-91           Zr-92           Zr-94
     > 89.9043D0,      90.9053D0,      91.9046D0,      93.9061D0,
*      Nb-93
     > 92.9060D0,      0.0,             0.0,             0.0,
*      Mo-98           Mo-92           Mo-95           Mo-96
     > 97.9055D0,      91.9063D0,      94.90584D0,     95.9046D0,
*      Tc
     > 98.0D0,         0.0,             0.0,             0.0,
*      Ru-102          Ru-99           Ru-100          Ru-104
     > 101.9037D0,     98.9061D0,      99.9030D0,      103.9055D0,
*      Rh-103
     > 102.9048D0,     0.0,             0.0,             0.0,
*      Pd-106          Pd-104           Pd-105         Pd-108
     > 105.9032D0,     103.9036D0,      104.9046D0,    107.90389D0,
*      Ag-107          Ag-109
     > 106.90509d0,    108.9047D0,      0.0,             0.0,
*      Cd-114          Cd-110           Cd-111         Cd-112
     > 113.9036D0,     109.9030D0,      110.9042D0,    111.9028D0,
*      In-115          In-113
     > 114.9041D0,     112.9043D0,      0.0,             0.0,
*      Sn-118          Sn-116           Sn-117         Sn-119
     > 117.9018D0,     115.9021D0,      116.9031D0,    118.9034D0,
*      Sb-121          Sb-123
     > 120.9038D0,     122.9041D0,      0.0,             0.0,
*      Te-130          Te-125           Te-126         Te-128
     > 129.9067D0,     124.9044D0,      125.9032D0,    127.9047D0,
*      I-127
     > 126.9004D0,     0.0,             0.0,             0.0,
*      Xe-132          Xe-129           Xe-131         Xe-134
     > 131.9042D0,     128.9048D0,      130.9051D0,    133.9054D0/
*
      end
*
