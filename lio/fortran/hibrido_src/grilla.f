      SUBROUTINE GRILLA(Q)
      use hibrido_common
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AR1(NPART),Q(NTQ)

C-----I2 ES EL Cl
      I2 = 1

      IF(3*NWAT.EQ.NPART) RETURN
c      write(*,*)' natom,nwat',natom,nwat,natom+nwat
C-----IDENTIFICO EL OXIGENO MAS CERCANO AL CL
      DO 12 I=3,NWAT+NATOM
c      write(*,*)'1grilla ',i,x(i),y(i),z(i)
c      write(*,*)'        ',i2,x(i2),y(i2),z(i2)
      IF(I.EQ.4.OR.I.EQ.5)GOTO 13
      AR1(I)=DSQRT((X(I)-X(I2))**2+(Y(I)-Y(I2))**2
     &         +(Z(I)-Z(I2))**2)
c      write(*,*)'2grilla r: ',i,i2,ar1(i)

      IF(I.EQ.3)THEN
      IOXIG=3
      GOTO 13
      ENDIF

      RZ=MIN(AR1(I),AR1(IOXIG))
      IF(RZ.EQ.AR1(I))THEN
      IOXIG=I
      ENDIF
13    CONTINUE
12    CONTINUE
c      write(*,*)'ioxig',ioxig

C-----EJE Z INSTANTANEO: O-CL , I2=CLORO
      XE=X(IOXIG)-X(I2)
      YE=Y(IOXIG)-Y(I2)
      ZE=Z(IOXIG)-Z(I2)
      EJEZ=DSQRT(XE**2+YE**2+ZE**2)
      XE=XE/EJEZ
      YE=YE/EJEZ
      ZE=ZE/EJEZ

C-----PROYECCION SOBRE EL EJE Z5 Y EL PERPENDICULAR R5
      DO 99 I=1,NPART
      Z5=(X(I)-X(I2))*XE+(Y(I)-Y(I2))*YE+(Z(I)-Z(I2))*ZE
      ZZ2=Z5*Z5
      RR2=(X(I)-X(I2))**2+(Y(I)-Y(I2))**2+(Z(I)-Z(I2))**2
      R5=DSQRT(RR2-ZZ2)
      IZZ = NINT(Z5/DELR2)
      IF(Z5.LT.ZERO)THEN
      IZZ = MAX(IZZ,-IGRIL/2)
      ELSE
      IZZ = MIN(IGRIL/2,IZZ)
      ENDIF
    
      IF(I.EQ.I2) IZZ = 0

      IZZ = IZZ + IGRIL/2 + 1

      IRR = NINT(R5/DELR2)
      IRR = MIN(IGRIL,IRR)
      IF(I.EQ.I2.OR.I.EQ.IOXIG) IRR=0
      IRR = IRR + 1
       IF(I.LE.NATOM)THEN
       QK(IZZ,IRR) = QK(IZZ,IRR) + Q(I)
       ELSE
       QK(IZZ,IRR) = QK(IZZ,IRR) + PC(I)/EE
       ENDIF
c       if(qk(izz,irr).ne.zero)write(*,*)'  qk '
c     &  ,izz,irr,qk(izz,irr)
C-----CL-O
      IF(I.EQ.3.OR.(I.GE.NATOM+1.AND.I.LE.NATOM+NWAT))THEN 
       M=1
       KG1(M,IZZ,IRR) = KG1(M,IZZ,IRR)+1
c       if(kg1(m,izz,irr).ne.0)write(*,*)' kg1 '
c     &  ,izz,irr,kg1(m,izz,irr)
C-----CL-EL RESTO (H)
      ELSE
       M=2
       KG1(M,IZZ,IRR) = KG1(M,IZZ,IRR)+1
c       if(kg1(m,izz,irr).ne.0)write(*,*)' kg1 '
c     &  ,izz,irr,kg1(m,izz,irr)
      ENDIF
      
99    CONTINUE 
      RETURN
      END
