      SUBROUTINE CVEMOTN1(NATSOL)
      INCLUDE 'COMM'
      DIMENSION AXX(NAT),AYY(NAT),AZZ(NAT)


      NOFSET = NSPECQ 
      DO 100 I = 1+NSPECQ, NSPECQ+NATSOL
      IF (I.GT.(1+NSPECQ)) NOFSET = NOFSET + NNAT(I-1)
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

      NOFSET = NSPECQ+1 
      DO 102 I = 1+NSPECQ, NSPECQ+NATSOL
      IF (I.GT.( 1+NSPECQ)) NOFSET = 
     > NOFSET + NNAT(I-1)
      DO 102 J= NOFSET,NNAT(I)+NOFSET 

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

C      NOFSET = 0
C      DO 1102 I = 1, NSPECQ+NATSOL
C      IF (I.GT.1) NOFSET = NOFSET + NNAT(I-1)
C      DO 1102 J= 1+NOFSET,NNAT(I)+NOFSET 

C      IF(ICLSTR.EQ.1)THEN 

C      IF(NSPECQ.EQ.0)THEN
C      XX(J) = XX(J) - XCM
C      YY(J) = YY(J) - YCM
C      ZZ(J) = ZZ(J) - ZCM
C      ENDIF
C      ENDIF
C      VX0(J) = VX0(J) - SX 
C      VY0(J) = VY0(J) - SY
C      VZ0(J) = VZ0(J) - SZ 

C      VX(J) = VX0(J)
C      VY(J) = VY0(J)
C      VZ(J) = VZ0(J)

C      CONTINUE

C      TEMPAV= TEMPAV + ( VX0(J)*VX0(J) + VY0(J)*VY0(J) + VZ0(J) *
C     &  VZ0(J) )*WWM(I)

c1102  CONTINUE

C-----TEMPERATURA DEL SISTEMA CUANTICO
c      NOFSET=0
c      DO 1202 I = 1, NSPECQ
c      IF (I.GT.1) NOFSET = NOFSET + NNAT(I-1)
c      DO 1202 J= 1+NOFSET,NNAT(I)+NOFSET 

c      TEMPSLT= TEMPSLT + ( VX0(J)*VX0(J) + VY0(J)*VY0(J) + VZ0(J) *
c     &  VZ0(J) )*WWM(I)

c1202  CONTINUE
       
C-----TEMPERATURA DEL SISTEMA CLASICO
      NOFSET=NATOM
      DO 1302 I = 1+NATOM , NATOM +NATSOL
      IF (I.GT.1+NATOM ) NOFSET = NOFSET + NNAT(I-1)
      DO 1302 J= 1+NOFSET,NNAT(I)+NOFSET 

      TEMPSLV= TEMPSLV + ( VX0(J)*VX0(J) + VY0(J)*VY0(J) + VZ0(J) *
     &  VZ0(J) )*WWM(I)

1302  CONTINUE
      TEMPAV=TEMPSLV

C=========temperatura constante

      tes= 1. + (DELTAT/KBER)*((TEMPRQ  /(TEMPSLV/GDFSLV)) -1.)
C     write (*,*) 'nanito', TEMPAV,TEMPAV/GDF,TEMPRQ,tes,GDF
      tes=sqrt(tes)
C-----UPDATE POSITIONS

      DO 12 I = 1+NSPECQ, NPART

      AX1 = X0(I)
      X0(I) = X(I)
      X(I) = tes*(XX(I)-X(I)) + X(I)
      XX(I) = AX1

      AY1= Y0(I)
      Y0(I) = Y(I)
      Y(I) = tes*(YY(I)-Y(I)) + Y(I)
      YY(I) = AY1

      AZ1= Z0(I)
      Z0(I) = Z(I)
      Z(I) = tes*(ZZ(I)-Z(I)) + Z(I)
      ZZ(I) = AZ1
 12   CONTINUE

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

