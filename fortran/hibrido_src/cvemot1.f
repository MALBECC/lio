      SUBROUTINE CVEMOT1(NATSOL)
      INCLUDE 'param'
      INCLUDE 'COMM'
      DIMENSION AXX(NAT),AYY(NAT),AZZ(NAT)
      
c       write(*,*)'entra a cve1'
c      ncn=0
c      write(*,*)'cve1'
c      do i=1+natom,nwat+natom
c      i1=i+nwat
c      i2=i1+nwat
c      d1=dsqrt((x(i)-x(i1))**2+(y(i)-y(i1))**2+(z(i)-z(i1))**2)
c      d2=dsqrt((x(i2)-x(i1))**2+(y(i2)-y(i1))**2+(z(i2)-z(i1))**2)
c      d3=dsqrt((x(i)-x(i2))**2+(y(i)-y(i2))**2+(z(i)-z(i2))**2)
c      if(d1+d2+d3-da(1)-da(2)-da(3).lt.1.D-06)ncn=ncn+1
c      enddo
c      if(ncn.eq.nwat)write(*,*)'OK'
c      if(ncn.ne.nwat)write(*,*)'NO-OK'

       NOFSET = 0
      DO 100 I = 1, NSPECQ+NATSOL
c        write(*,*)'fff ',i,fff(i)
      IF (I.GT.1) NOFSET = NOFSET + NNAT(I-1)
      DO 100 J= 1+NOFSET,NNAT(I)+NOFSET 
c--- Si no hizo un calculo DFT:
        if(j.le.natom)then
        xx(j) = x(j)
        yy(j) = y(J)
        zz(j) = z(j)
        else
        XX(J) = TWO * X(J) -  X0(J) + FFF(I) * FX(J)
        YY(J) = TWO * Y(J) -  Y0(J) + FFF(I) * FY(J)
        ZZ(J) = TWO * Z(J) -  Z0(J) + FFF(I) * FZ(J)
        endif
      AXX(J) = XX(J)
      AYY(J) = YY(J)
      AZZ(J) = ZZ(J)
c      write(*,*)'1-cve1 xx fff fx',j,xx(j),fff(i),fx(j)
100   CONTINUE
 
c      ncn=0
c      write(*,*)'cve2'
c      do i=1+natom,nwat+natom
c      i1=i+nwat
c      i2=i1+nwat
c      d1=dsqrt((xx(i)-xx(i1))**2+(yy(i)-yy(i1))**2+(zz(i)-zz(i1))**2)
c      d2=dsqrt((xx(i2)-xx(i1))**2+(yy(i2)-yy(i1))**2+(zz(i2)-zz(i1))**2)
c      d3=dsqrt((xx(i)-xx(i2))**2+(yy(i)-yy(i2))**2+(zz(i)-zz(i2))**2)
c      if(d1+d2+d3-da(1)-da(2)-da(3).lt.1.D-06)ncn=ncn+1
c      enddo
c      if(ncn.eq.nwat)write(*,*)'OK'
c      if(ncn.ne.nwat)write(*,*)'NO-OK'

      CALL GAMMA(NATSOL) 
      do i=1+natom,npart
c      write(*,*)'2-cve1 xx fff fx',i,xx(i)
      enddo
c      ncn=0
c      write(*,*)'cve3'
c      do i=1+natom,nwat+natom
c      i1=i+nwat
c      i2=i1+nwat
c      d1=dsqrt((xx(i)-xx(i1))**2+(yy(i)-yy(i1))**2+(zz(i)-zz(i1))**2)
c      d2=dsqrt((xx(i2)-xx(i1))**2+(yy(i2)-yy(i1))**2+(zz(i2)-zz(i1))**2)
c      d3=dsqrt((xx(i)-xx(i2))**2+(yy(i)-yy(i2))**2+(zz(i)-zz(i2))**2)
c      if(d1+d2+d3-da(1)-da(2)-da(3).lt.1.D-06)ncn=ncn+1
c      enddo
c      if(ncn.eq.nwat)write(*,*)'OK'
c      if(ncn.ne.nwat)write(*,*)'NO-OK'
 
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
     
      if(j.le.natom)then
      VX0(J) = ZERO 
      VY0(J) = ZERO 
      VZ0(J) = ZERO

      else
     
      VX0(J) = FACTV* (XX(J) - X0(J))
      VY0(J) = FACTV* (YY(J) - Y0(J))
      VZ0(J) = FACTV* (ZZ(J) - Z0(J))
       
      endif

c      write(*,*)'1-cve1 xx x0',j,xx(j),yy(j),zz(j),x0(j),y0(j),z0(j)
c      XCM = XCM + XX(J)*WWM(I)
c      YCM = YCM + YY(J)*WWM(I)
c      ZCM = ZCM + ZZ(J)*WWM(I)

c      SX = SX + VX0(J)*WWM(I)
c      SY = SY + VY0(J)*WWM(I)
c      SZ = SZ + VZ0(J)*WWM(I)

102   CONTINUE

c      XCM = XCM /TOTMAS
c      YCM = YCM /TOTMAS
c      ZCM = ZCM /TOTMAS

c      SX = SX /TOTMAS
c      SY = SY /TOTMAS
c      SZ = SZ /TOTMAS


C-----TEMPERATURA DEL SISTEMA ENTERO

      NOFSET = 0
      DO 1102 I = 1, NSPECQ+NATSOL
      IF (I.GT.1) NOFSET = NOFSET + NNAT(I-1)
      DO 1102 J= 1+NOFSET,NNAT(I)+NOFSET 

c        IF(NDFT.NE.1)THEN
c      XX(J) = XX(J) - XCM
c      YY(J) = YY(J) - YCM
c      ZZ(J) = ZZ(J) - ZCM
c      VX0(J) = VX0(J) - SX 
c      VY0(J) = VY0(J) - SY
c      VZ0(J) = VZ0(J) - SZ 
c        ELSE

      VX(J) = VX0(J)
      VY(J) = VY0(J)
      VZ(J) = VZ0(J)
c      write(*,*)'2-cve1 vx0',j,vx0(j),vy0(j),vz0(j)
c      pause

      CONTINUE
c        ENDIF


      TEMPAV= TEMPAV + ( VX0(J)*VX0(J) + VY0(J)*VY0(J) + VZ0(J) *
     &  VZ0(J) )*WWM(I)
c      write(*,*)'cve1 tempav',j,tempav

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

c      write(*,*)'sale de cve1'

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CVTMOT1(NATSOL)
      INCLUDE 'param'
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
      if(j.le.natom)then
      xx(j) = x(j)
      yy(j) = y(j)
      zz(j) = z(j)
      endif
      XX(J) = TWO * X(J) -  X0(J) + FFF(II) * FX(J) -SLD * VX(J)
      YY(J) = TWO * Y(J) -  Y0(J) + FFF(II) * FY(J) -SLD * VY(J)
      ZZ(J) = TWO * Z(J) -  Z0(J) + FFF(II) * FZ(J) -SLD * VZ(J)
4100  CONTINUE
10    CONTINUE

c         IF(NDFT.EQ.1)THEN
c       DO 77 I=1,NATOM
c       X0(I)=X(I)
c       Y0(I)=Y(I)
c       Z0(I)=Z(I)
c       XX(I)=X(I)
c       YY(I)=Y(I)
c       ZZ(I)=Z(I)
c77     CONTINUE       
c         ELSE
c        CONTINUE       
c         ENDIF

      CALL GAMMA(NATSOL)

      DO 11 I = 1, NPART 
      if(i.le.natom)then
      vx0(i)=zero
      vy0(i)=zero
      vz0(i)=zero
      else 
      
      VX0(I) = FACTV*(XX(I) - X0(I))
      VY0(I) = FACTV*(YY(I) - Y0(I))
      VZ0(I) = FACTV*(ZZ(I) - Z0(I))
      endif
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

        IF(NDFT.EQ.1)THEN
      DO 229 I = 1, NPART
      VX0(I) = FACTV*(XX(I) - X0(I))  
      VY0(I) = FACTV*(YY(I) - Y0(I))
      VZ0(I) = FACTV*(ZZ(I) - Z0(I)) 
229   CONTINUE

        ELSE

      DO 228 I = 1, NPART
      XX(I) = XX(I) - XXCM
      YY(I) = YY(I) - YYCM
      ZZ(I) = ZZ(I) - ZZCM

      VX0(I) = FACTV*(XX(I) - X0(I))  
      VY0(I) = FACTV*(YY(I) - Y0(I))
      VZ0(I) = FACTV*(ZZ(I) - Z0(I)) 

228   CONTINUE

        ENDIF

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


