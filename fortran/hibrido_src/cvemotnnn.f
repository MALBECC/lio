      SUBROUTINE CVEMOTNNN(NATSOL)

C Esta subrutina mueve los atomos del soluto, dejando fijos 
C los atomos de solvente.

      INCLUDE 'param'
        INCLUDE 'COMM'
      DIMENSION AXX(NAT),AYY(NAT),AZZ(NAT)

      NOFSET = 0
      DO 100 I = 1, NSPECQ+NATSOL
      IF (I.GT.1) NOFSET = NOFSET + NNAT(I-1)
      DO 100 J= 1+NOFSET,NNAT(I)+NOFSET 
      XX(J) = TWO * X(J) -  X0(J) + FFF(I) * FX(J)
      YY(J) = TWO * Y(J) -  Y0(J) + FFF(I) * FY(J)
      ZZ(J) = TWO * Z(J) -  Z0(J) + FFF(I) * FZ(J)
      AXX(J) = XX(J)
      AYY(J) = YY(J)
      AZZ(J) = ZZ(J)
100   CONTINUE
 
      CALL GAMMA(NATSOL)
      SX = ZERO
      SY = ZERO
      SZ = ZERO

      TEMPAV  = ZERO
      TEMPSLV = ZERO
      TEMPSLT = ZERO

      XCM =ZERO
      YCM =ZERO
      ZCM =ZERO

      NOFSET = 0
      DO 102 I = 1, NSPECQ+NATSOL
      IF (I.GT.1) NOFSET = NOFSET + NNAT(I-1)
      DO 102 J= 1+NOFSET,NNAT(I)+NOFSET 

      VX0(J) = FACTV* (XX(J) - X0(J))
      VY0(J) = FACTV* (YY(J) - Y0(J))
      VZ0(J) = FACTV* (ZZ(J) - Z0(J))
      
      XCM = XCM + XX(J)*WWM(I)
      YCM = YCM + YY(J)*WWM(I)
      ZCM = ZCM + ZZ(J)*WWM(I)

      SX = SX + VX0(J)*WWM(I)
      SY = SY + VY0(J)*WWM(I)
      SZ = SZ + VZ0(J)*WWM(I)

102   CONTINUE

      XCM = XCM /TOTMAS
      YCM = YCM /TOTMAS
      ZCM = ZCM /TOTMAS

      SX = SX /TOTMAS
      SY = SY /TOTMAS
      SZ = SZ /TOTMAS


C-----TEMPERATURA DEL SISTEMA ENTERO

      NOFSET = 0
      DO 1102 I = 1, NSPECQ+NATSOL
      IF (I.GT.1) NOFSET = NOFSET + NNAT(I-1)
      DO 1102 J= 1+NOFSET,NNAT(I)+NOFSET 

      IF(ICLSTR.EQ.1)THEN 

      IF(NSPECQ.EQ.0)THEN
      XX(J) = XX(J) - XCM
      YY(J) = YY(J) - YCM
      ZZ(J) = ZZ(J) - ZCM
      ENDIF
      ENDIF
      VX0(J) = VX0(J) - SX 
      VY0(J) = VY0(J) - SY
      VZ0(J) = VZ0(J) - SZ 

      VX(J) = VX0(J)
      VY(J) = VY0(J)
      VZ(J) = VZ0(J)

      CONTINUE

      TEMPAV= TEMPAV + ( VX0(J)*VX0(J) + VY0(J)*VY0(J) + VZ0(J) *
     &  VZ0(J) )*WWM(I)

1102  CONTINUE

C-----TEMPERATURA DEL SISTEMA CUANTICO
      NOFSET=0
      DO 1202 I = 1, NSPECQ
      IF (I.GT.1) NOFSET = NOFSET + NNAT(I-1)
      DO 1202 J= 1+NOFSET,NNAT(I)+NOFSET 

      TEMPSLT= TEMPSLT + ( VX0(J)*VX0(J) + VY0(J)*VY0(J) + VZ0(J) *
     &  VZ0(J) )*WWM(I)

1202  CONTINUE
       
C-----TEMPERATURA DEL SISTEMA CLASICO
      NOFSET=NATOM
      DO 1302 I = 1+NATOM , NATOM +NATSOL
      IF (I.GT.1+NATOM ) NOFSET = NOFSET + NNAT(I-1)
      DO 1302 J= 1+NOFSET,NNAT(I)+NOFSET 

      TEMPSLV= TEMPSLV + ( VX0(J)*VX0(J) + VY0(J)*VY0(J) + VZ0(J) *
     &  VZ0(J) )*WWM(I)

1302  CONTINUE
C=========temperatura constante
C=========Solvente
      tes= 1. + (DELTAT/kber)*((TEMPRSV  /(TEMPSLV/GDFSLV)) -1.)
C      write (*,*) 'nanito', TEMPSLV/GDFSLV,TEMPRQ,tes
      tes=sqrt(tes)
C-----UPDATE POSITIONS

C "NO HACE UPDATE DE LAS POSICIONES CLASICAS"

      DO 12 I = NSPECQ + 1, NPART
 
      AX1 = X0(I)
      X0(I) = X(I)
c     X(I) = tes*(XX(I)-X(I)) + X(I)
      XX(I) = AX1
c
      AY1= Y0(I)
      Y0(I) = Y(I)
c     Y(I) = tes*(YY(I)-Y(I)) + Y(I)
      YY(I) = AY1
c
      AZ1= Z0(I)
      Z0(I) = Z(I)
c     Z(I) = tes*(ZZ(I)-Z(I)) + Z(I)
      ZZ(I) = AZ1
 12   CONTINUE

C=========Soluto  
c      write(*,*) 'KBER=',KBER
      tess= 1. + (DELTAT/KBER)*((TEMPRQ  /(TEMPSLT/GDFSLT)) -1.)
c      write (*,*) 'nanito',TEMPRQ,DELTAT,
c     >TEMPSLT/GDFSLT,tess
      tess=sqrt(tess)
C-----UPDATE POSITIONS
 
      DO 122 I = 1, NSPECQ 
 
      AX1 = X0(I)
      X0(I) = X(I)
      X(I) = tess*(XX(I)-X(I)) + X(I)
      XX(I) = AX1
 
      AY1= Y0(I)
      Y0(I) = Y(I)
      Y(I) = tess*(YY(I)-Y(I)) + Y(I)
      YY(I) = AY1
 
      AZ1= Z0(I)
      Z0(I) = Z(I)
      Z(I) = tess*(ZZ(I)-Z(I)) + Z(I)
      ZZ(I) = AZ1
 122   CONTINUE

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

