      SUBROUTINE CVEMOT(NATSOL)
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


C-----UPDATE POSITIONS

      DO 12 I = 1, NPART

      AX1 = X0(I)
      X0(I) = X(I)
      X(I) = XX(I)
      XX(I) = AX1

      AY1= Y0(I)
      Y0(I) = Y(I)
      Y(I) = YY(I)
      YY(I) = AY1

      AZ1= Z0(I)
      Z0(I) = Z(I)
      Z(I) = ZZ(I)
      ZZ(I) = AZ1
 12   CONTINUE

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CVTMOT(NATSOL)
      INCLUDE 'COMM'


      DO 2922 ITER = 1, 200

      TEMPAV  = ZERO
      TEMPSLV = ZERO
      TEMPSLT = ZERO

      NOFSET = 0

      DO 145 II =1 ,NSPECQ+NATSOL
      IF (II.GT.1) NOFSET = NOFSET + NNAT(II-1)
      DO 145 I= 1+NOFSET, NNAT(II)+NOFSET
      TEMPAV = TEMPAV + (VX(I)*VX(I) + VY(I)*VY(I) + VZ(I)*VZ(I))
     & * WWM(II)
 145  CONTINUE

C-----DYNAMICS OF TIME SCALING

      SN = TWO * ST - S0 + DEL/QS*(TEMPAV - GDF*TEMPRQ)
      SD0 = FACTV * ( SN - S0)
      SLD = SD0*DEL
     

      NOFSET = 0
      DO 10 II =1 ,NSPECQ+NATSOL
      IF (II.GT.1) NOFSET = NOFSET + NNAT(II-1)
      DO 4100 J= 1+NOFSET, NNAT(II)+NOFSET
      XX(J) = TWO * X(J) -  X0(J) + FFF(II) * FX(J) -SLD * VX(J)
      YY(J) = TWO * Y(J) -  Y0(J) + FFF(II) * FY(J) -SLD * VY(J)
      ZZ(J) = TWO * Z(J) -  Z0(J) + FFF(II) * FZ(J) -SLD * VZ(J)
4100  CONTINUE
10    CONTINUE


      CALL GAMMA(NATSOL)

      DO 11 I = 1, NPART 
      
      VX0(I) = FACTV*(XX(I) - X0(I))
      VY0(I) = FACTV*(YY(I) - Y0(I))
      VZ0(I) = FACTV*(ZZ(I) - Z0(I))
11    CONTINUE

      DO 111 I = 1, NPART
      IF (VX0(I).NE.ZERO)THEN
      DIFF = DABS ( (VX0(I) - VX(I)) / VX0(I) )
      IF (DIFF.GT.0.0001)GO TO 383
      ENDIF
      IF (VY0(I).NE.ZERO)THEN
      DIFF = DABS ( (VY0(I) - VY(I)) / VY0(I) )
      IF (DIFF.GT.0.0001)GO TO 383
      ENDIF 
      IF (VZ0(I).NE.ZERO)THEN
      DIFF = DABS ( (VZ0(I) - VZ(I)) / VZ0(I) )
      IF (DIFF.GT.0.0001)GO TO 383
      ENDIF 
111   CONTINUE

      GOTO 66

383   CONTINUE

      DO 2881 I = 1, NPART
      VX(I) = VX0(I)
      VY(I) = VY0(I)
      VZ(I) = VZ0(I)
2881  CONTINUE



2922  CONTINUE

      WRITE (6,*) ' TOO MANY ITERATIONS IN CVTMOT'


      STOP
     
C-----NEW ESTIMATION FOR THE NEXT STEP
66    NOFSET = 0
      FACT2 = ONE/DELTAT
      DO 33 I = 1, NSPECQ+NATSOL
      IF (I.GT.1)NOFSET = NOFSET + NNAT(I-1)
      FACTC = FFF(I) / DELTAT 
      DO 533 J = 1+NOFSET,NNAT(I)+NOFSET
      VX(J) = VX0(J) + FACTC*FX(J) - FACT2*(SLD*VX(J)-GX(J))
      VY(J) = VY0(J) + FACTC*FY(J) - FACT2*(SLD*VY(J)-GY(J))
      VZ(J) = VZ0(J) + FACTC*FZ(J) - FACT2*(SLD*VZ(J)-GZ(J))
      
533   CONTINUE
33    CONTINUE

      XXCM =ZERO
      YYCM =ZERO
      ZZCM =ZERO

      IOFSET = 0
      DO 102 I = 1, NSPECQ+NATSOL
      IF (I.GT.1) IOFSET = IOFSET + NNAT(I-1)
      DO 8101 J= 1+IOFSET, IOFSET+NNAT(I)

      XXCM = XXCM + XX(J)*WWM(I)
      YYCM = YYCM + YY(J)*WWM(I)
      ZZCM = ZZCM + ZZ(J)*WWM(I)

8101  CONTINUE
102   CONTINUE

      XXCM = XXCM /TOTMAS
      YYCM = YYCM /TOTMAS
      ZZCM = ZZCM /TOTMAS


      DO 228 I = 1, NPART
c      XX(I) = XX(I) - XXCM
c      YY(I) = YY(I) - YYCM
c      ZZ(I) = ZZ(I) - ZZCM

      VX0(I) = FACTV*(XX(I) - X0(I))  
      VY0(I) = FACTV*(YY(I) - Y0(I))
      VZ0(I) = FACTV*(ZZ(I) - Z0(I)) 

228   CONTINUE


      IOFSET = 0
      TEMPAV = ZERO

      DO 8288 I=1,NSPECQ+NATSOL
      TEMPA(I) = ZERO
8288  CONTINUE

C-----TEMPERATURA DEL SISTEMA ENTERO      
      DO 1102 I = 1, NSPECQ+NATSOL
      IF (I.GT.1) IOFSET = IOFSET + NNAT(I-1)
      DO 1101 J = IOFSET+1, IOFSET+NNAT(I)
      TEMPA(I)= TEMPA(I)+(VX0(J)*VX0(J)+VY0(J)*VY0(J)+VZ0(J)*
     &   VZ0(J))*WWM(I)
1101  CONTINUE
      TEMPAV = TEMPAV + TEMPA(I)
1102  CONTINUE

C-----TEMPERATURA DEL SOLUTO             
      NOFSET=0
      DO 1202 I = 1, NSPECQ
      IF (I.GT.1) NOFSET = NOFSET + NNAT(I-1)
      DO 1202 J= 1+NOFSET,NNAT(I)+NOFSET 

      TEMPSLT= TEMPSLT + ( VX0(J)*VX0(J) + VY0(J)*VY0(J) + VZ0(J) *
     &  VZ0(J) )*WWM(I)

1202  CONTINUE
       
C-----TEMPERATURA DEL SOLVENTE        
      NOFSET=NATOM
      DO 1302 I = 1+NSPECQ, NSPECQ+NATSOL
      IF (I.GT.1+NSPECQ) NOFSET = NOFSET + NNAT(I-1)
      DO 1302 J= 1+NOFSET,NNAT(I)+NOFSET 

      TEMPSLV= TEMPSLV + ( VX0(J)*VX0(J) + VY0(J)*VY0(J) + VZ0(J) *
     &  VZ0(J) )*WWM(I)

1302  CONTINUE




C-----UPDATE POSITIONS

      DO 12 I = 1, NPART
      X2 = X0(I)
      X0(I) = X(I)
      X(I) = XX(I)
      XX(I) = X2 

      Y2= Y0(I)
      Y0(I) = Y(I)
      Y(I) = YY(I)
      YY(I) = Y2

      Z2= Z0(I)
      Z0(I) = Z(I)
      Z(I) = ZZ(I)
      ZZ(I) = Z2

 12   CONTINUE


      S0 = ST
      ST = SN

      RETURN
      END


