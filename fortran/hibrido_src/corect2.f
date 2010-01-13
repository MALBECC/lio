      SUBROUTINE CORECT2(NATSOL,RM,EM,NTQ,NSS,PM,NT,AVNPU,
     & PMAX,PZMAX,HISTO,HISTO2,DELP,DELPZ)
      INCLUDE 'COMM'
      DIMENSION EM(NTQ+NSS),RM(NTQ+NSS),PM(NT),HISTO(100,100),
     & HISTO2(100)


      E2 = EE*EE
C---------------------------------------------------------------------*
C                    (1) 00, (2)OH , (3)HH
C                                                         
C---------------------------------------------------------------------*
c     write(*,*) 'Estoy en corect2'
C-----PASO A UNIDADES PROGRAMA LOS PARAM L-J CUANTICOS
      DO I=1,NSPECQ
      IF(EM(I).LT.ZERO) EM(I)=-EM(I)
      EPS(I) = EM(I)*(4.3598D05)/BOLTZF
      SIGMA(I) = RM(I)*0.52918D0
c     WRITE (*,*) 'SIG=',SIGMA(I)
c     WRITE (*,*) 'EPS=',EPS(I)
      ENDDO
      
      AZ(1) =  ZZZ(1)*ZZZ(1)*E2
      AZ(2) =  ZZZ(1)*ZZZ(2)*E2
      AZ(3) =  ZZZ(2)*ZZZ(2)*E2

C-----LJ CLASICOS
      DO I=1+NSPECQ,NSPECQ+NATSOL
      EPSI = EPS(I) * FOUR * BOLTZF
      SIG3 = SIGMA(I)
      SIG3 = SIG3 * SIG3 * SIG3
      SIG6 = SIG3 * SIG3
      SIG12 = SIG6 * SIG6
      E6(I) = EPSI * SIG6
      E12(I) = EPSI * SIG12
      F12(I) = 12.D0 * E12(I)
      F6(I) = 6.D0 * E6(I)
      
      ENDDO

C-----PARAM LENN-JONES OXIGENOS CLASICOS
      EE12O=E12(NSPECQ+1)
      EE6O=E6(NSPECQ+1)
      FF12O=F12(NSPECQ+1)
      FF6O=F6(NSPECQ+1)

C-----PARAM LENN-JONES HIDOGENOS CLASICOS
      EE12H=E12(NSPECQ+2)
      EE6H=E6(NSPECQ+2)
      FF12H=F12(NSPECQ+2)
      FF6H=F6(NSPECQ+2)

C-----PARAM LENN-JONES HID-OX
        SIG3OH  = PTFIVE * (SIGMA(NSPECQ+2)+SIGMA(NSPECQ+1))
        EPSSOH  = DSQRT (EPS(NSPECQ+2)*EPS(NSPECQ+1))
        IF (SIGMA(I).EQ.ZERO) THEN
        SIG3OH = ZERO
        EPSSOH = ZERO
        ENDIF

        EPSIOH  = EPSSOH * FOUR * BOLTZF
        SIG3OH  = SIG3OH * SIG3OH * SIG3OH
        SIG6OH  = SIG3OH * SIG3OH
        SIG12OH = SIG6OH * SIG6OH
        EE6OH   = EPSI * SIG6OH
        EE12OH  = EPSI * SIG12OH
        FF12OH  = 12.D0 * EE12OH
        FF6OH   = 6.D0  * EE6OH



c      IF (IEWLD.NE.1)THEN
c      RCT = DBLE((BXLGTH/TWO)*XFAX)
c      RCU1 = RCT - RBUF
c      RM2 = RCU1**2
c      write(*,*) 'aca va RM2',RM2,'=',RCU1,'**2',RCT,RBUF
c      UC = ONE/RCU1
c      FC = ONE / RM2
c      B = FC / (4.D0*RCU1*(RM2-RCTSQ))
c      AB = -2.D0*B*RCTSQ
c      CB = UC + (B*RM2+AB)*RM2
c      CSHIFT = (AB+B*RCTSQ)*RCTSQ
c      ESHIFT = -CB + CSHIFT
c      AF = TWO*AB
c      BF = 4.D0*B
c      write(*,*) 'aca va CSHIFT y AB',CSHIFT,AB
c      ENDIF
C--------------------------------------------------C
C       [OV]=(erg-16/e2)
C       [EN]=(erg-16/e)
C       OV(1)= OO (R=0)      OV(3)= HH (R=0)
C       OV(2)= OH            OV(4)= HH    
C       ENOH= Electronegativ(O) - Electronegativ(H)
C---------------------------------------------------C       
      ENOH =47581.715D0 /EE
      OV(1)=258159.81D0/E2
      OV(2)=198969.24D0/E2
      OV(3)=245237.93D0/E2
      OV(4)=141446.01D0/E2
      DENOM=TWO*OV(1)+OV(3)-FOUR*OV(2)+OV(4)
      EGP=-ENOH*ENOH/DENOM
      AA1=ENOH
      BB=(OV(1)+OV(3))*PTFIVE-OV(2)
      DD=OV(3)-OV(4)


C-----CALCULA COSAS QUE NECESITA DESPUES
      IDUM = 0
      NATOM = 0
      DO I=1,NSPECQ
      NATOM=NATOM+NNAT(I)
      ENDDO
      NWAT = NNAT(NSPECQ+1)
      NVINC = 3
      NPART = 0
      DO 121 I = NSPECQ+1, NSPECQ+NATSOL
      NPART = NPART + NNAT(I)
121   CONTINUE
      NPART = NPART + NATOM
      NDSLV = NWAT*(9 - NVINC)
      NDGREE = NDSLV
      GDFSLV= DBLE(NDGREE)*BOLTZF
c       write(*,*)'GDFSLT',GDFSLT
      IF(NDFT.NE.1)THEN
        GDFSLT = ZER0
        GDFSLT = 0.0D0
        caca = ZERO
c       write(*,*) 'CACA',CACA
c       write(*,*)'GDFSLT',GDFSLT
      ELSE
        GDFSLT= DBLE(3*NATOM)*BOLTZF
      ENDIF
      GDF = GDFSLV + GDFSLT - 6.D0*BOLTZF
c     write(*,*) 'nanitocorect',GDF,GDFSLT,GDFSLV
      GDQ = DBLE(2*NWAT) * BOLTZF
      VFTR = DSQRT(BOLTZF*TEMPRQ)
      RCT = DBLE((BXLGTH/TWO)*XFAX)
      RCTSQ = RCT * RCT
      RKAPPA =  DBLE(5.D0/BXLGTH)
      RKAPPA2 = RKAPPA*RKAPPA
      QB = ONE/DBLE(4*RKAPPA2)
      QA = TWOPI/(BXLGTH**3)
      KSQMAX = KMAX*KMAX
      AXI = ONE/BXLGTH
      BYI = ONE/BXLGTH
      CZI = ONE/BXLGTH
      ECKT = ONE / (TEMPRQ*BOLTZF)
      CC = ECKT / DBLE(NPART)
      QS = QS * GDF
      ST = ZERO
      S0 = ZERO
      SD = ZERO
      SS = ZERO
      TOLL= 1.D-07
      CONVF=ONE/694.725D0

      ECKTQ = ONE / (TEMPRQQ*BOLTZF)
      CCQ = ECKTQ / DBLE(NPART)
      QSQ = QSQ * GDQ
      STQ = ZERO
      S0Q = ZERO
      SDQ = ZERO
      SSQ = ZERO

      MAXIT = 400
      MAXI5 = 200
      TEMPS = ZERO
      AVNPU = AVNUM * 1.D-24
      FACTA = ONE/(THREE*BOLTZF*DBLE(NPART))
      FACTV = ONE / (TWO*DELTAT)
      DEL = DELTAT * DELTAT
      REDEL = ONE / DELTAT
      FACT2 =  TWO / DELTAT
      WEIGHT = ZERO
      TOTMAS = ZERO
      SUMFX = ZERO
      SUMFY = ZERO
      SUMFZ = ZERO
      SUMF  = ZERO
      ENEFUR = ZERO

C-----TOTMAS=MASA TOTAL 
      DO I = 1, NSPECQ+NATSOL
      WWM(I) = WWM(I) / AVNPU
      WEIGHT = WEIGHT + WWM(I)
      TOTMAS = TOTMAS + WWM(I) * NNAT(I)
      VF(I) = VFTR / DSQRT(WWM(I))
      FFF(I) = DEL / WWM(I)
      ENDDO
      WWQ = WWM(NSPECQ+NATSOL+1)
      TOTMASQ = DBLE(NWAT)*THREE*WWQ
      VFTRQ = DSQRT(BOLTZF*TEMPRQQ)
      VFQ =  VFTRQ / DSQRT(WWQ)
      FFQ =  DEL / WWQ

C-----MASAS DE NATOM PARA INTSG 
C-----(NO INTERESAN LAS UNID. SINO LA RELACION PM(I)/MTOT)
      DO I=1,NATOM
      PM(I) = WWM(I)
      ENDDO

C-----MI NUMERACION PARA EL SOLVENTE EMPIEZA EN NATOM
      IF(ICON.EQ.0)THEN
      NOFSET=NATOM
      DO 55 I=NSPECQ+1,NSPECQ+NATSOL
      IF(I.GT.NSPECQ+1)NOFSET=NOFSET+NNAT(I-1)
      DO 55 J=NOFSET+1,NNAT(I)+NOFSET
      PC(J)=ZZZ(I-NSPECQ)*EE
      PC0(J)=PC(J)
55    CONTINUE
      ENDIF

      DO I=1,NATOM
      PC1(I)=PC1(I)*EE
      ENDDO

      DO 180 I = 1, 500
180   AC(I) = ZERO

      DO 10 I = 1, 8
      DO 10 J = 1, IGGRID
      KG(I,J) = 0
10    CONTINUE

      DO 33 I = 1 ,IGRIL+1
      DO 33 J = 1 ,IGRIL+1
      QK(I,J) = ZERO
      DO 33 M = 1,2
      KG1(M,I,J) = 0
33    CONTINUE

      DO 36 I=1,100
      DO 35 J=1,100
      HISTO(J,I)=0.0
35    CONTINUE
      HISTO2(I)=0.0
36    CONTINUE

      DELP = PMAX/100.d0
      DELPZ= PZMAX/100.d0

C-----ADIMENSIONAL
      IF(IUNID.EQ.1) FACTOR=CC
      
C-----EN KELVIN
      IF(IUNID.EQ.2) FACTOR=ONE/BOLTZF
     
C-----EN Kcal/mol
      IF(IUNID.EQ.3) FACTOR=CONVF

C-----EN HARTREES
      IF(IUNID.EQ.4) FACTOR=0.52918D0/E2

        DO I=NATOM+1,NATOM+NWAT
        AT(I)='O    '
        AT(I+NWAT)  ='H    '
        AT(I+2*NWAT)='H    '
        ENDDO

c======esto lo puso nano=== 
      IF (IEWLD.NE.1)THEN
      RCT = DBLE((BXLGTH/TWO)*XFAX)
      RCU1 = RCT - RBUF
      RM2 = RCU1**2
c      write(*,*) 'aca va RM2',RM2,'=',RCU1,'**2',RCT,RBUF
      UC = ONE/RCU1
      FC = ONE / RM2
      B = FC / (4.D0*RCU1*(RM2-RCTSQ))
      AB = -2.D0*B*RCTSQ
      CB = UC + (B*RM2+AB)*RM2
      CSHIFT = (AB+B*RCTSQ)*RCTSQ
      ESHIFT = -CB + CSHIFT
      AF = TWO*AB
      BF = 4.D0*B
C      write(*,*) 'aca va CSHIFT y AB',CSHIFT,AB
c==========HASTA ACA
      ENDIF

      RETURN
      END

