      SUBROUTINE GAMMA(NATSOL)
      INCLUDE 'param'
        INCLUDE 'COMM'
      LOGICAL FLAG
      DIMENSION DXT(NAT),DYT(NAT),DZT(NAT)
      DIMENSION DXP(NAT),DYP(NAT),DZP(NAT)
      DIMENSION DP2(NAT),WWMM(NPART)
      
c      write(*,*)' wwm ini ',wwm(nspecq+1),wwm(nspecq+2),wwm(nspecq+3)
      WWMM(1)=WWM(NSPECQ+1)
      WWMM(2)=WWM(NSPECQ+2)
      WWMM(3)=WWM(NSPECQ+3)
 

      if(1.ne.2)goto 8898
      write(*,*)'gamma',da(1),da(2),da(3),toll,maxi5,maxit,natom
      do i=1+natom,nwat+natom
      i1=i+nwat
      i2=i1+nwat
      d1=dsqrt((x(i)-x(i1))**2+(y(i)-y(i1))**2+(z(i)-z(i1))**2)
      d2=dsqrt((x(i2)-x(i1))**2+(y(i2)-y(i1))**2+(z(i2)-z(i1))**2)
      d3=dsqrt((x(i)-x(i2))**2+(y(i)-y(i2))**2+(z(i)-z(i2))**2)
      dt1=dsqrt((xx(i)-xx(i1))**2+(yy(i)-yy(i1))**2+
     & (zz(i)-zz(i1))**2)
      dt2=dsqrt((xx(i2)-xx(i1))**2+(yy(i2)-yy(i1))**2+
     & (zz(i2)-zz(i1))**2)
      dt3=dsqrt((xx(i)-xx(i2))**2+(yy(i)-yy(i2))**2+
     & (zz(i)-zz(i2))**2)
      write(*,*)'entra a gamma' ,i,d1,d2,d3,dt1,dt2,dt3
      write(*,*)'entra g',abs(dt1-d1),abs(dt2-d2),abs(dt3-d3)
      enddo
8898   continue
c      write(40,*)npart
c      write(40,*)
      DO 2221 I = NATOM+1, NPART
      GX(I) = XX(I)
      GY(I) = YY(I)
      GZ(I) = ZZ(I)
c      write(40,110)at(i),x(i),y(i),z(i),i
2221  CONTINUE
110   format(1x,a5,3F10.5,2x,i5)

      IF(NWAT.EQ.0)RETURN

      ITER = 0

889   NOFSET = NATOM
      DO 11 NN=1,NWAT
      NN1= NN + NOFSET
      NN2=NN1+NWAT
      NN3=NN2+NWAT
      N1 = NN + NWAT
      N2 = N1 + NWAT
      DXT(NN)=X(NN1)-X(NN2)
      DYT(NN)=Y(NN1)-Y(NN2)
      DZT(NN)=Z(NN1)-Z(NN2)
      DXT(N1)=X(NN1)-X(NN3)
      DYT(N1)=Y(NN1)-Y(NN3)
      DZT(N1)=Z(NN1)-Z(NN3)
      DXT(N2)=X(NN2)-X(NN3)
      DYT(N2)=Y(NN2)-Y(NN3)
      DZT(N2)=Z(NN2)-Z(NN3)
11    CONTINUE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     BEGINNING OF THE LOOP      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FACTO1=TWO*WWMM(1)*WWMM(2)/(FOUR*(WWMM(1)+WWMM(2)))
      FACTO2=TWO*WWMM(1)*WWMM(2)/(FOUR*(WWMM(1)+WWMM(2)))
      FACTO3=TWO*WWMM(2)*WWMM(2)/(FOUR*(WWMM(2)+WWMM(2)))

c       write(*,*)'facto 1',facto1,facto2,facto3
c       write(*,*)'wwmm 1',wwmm(1),wwmm(2),wwmm(3),two,four

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     BOND CONSTRAINTS           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

8000  FLAG=.TRUE.
      AD2M=ZERO
      
      NOFSET=NATOM
      DO 12 NN = 1, NWAT
      NN1= NOFSET + NN
      NN2= NN1+NWAT
      NN3= NN2+NWAT
      N1 = NN + NWAT
      N2 = N1 + NWAT
      DXP(NN)=XX(NN1)-XX(NN2)
      DYP(NN)=YY(NN1)-YY(NN2)
      DZP(NN)=ZZ(NN1)-ZZ(NN2)
      DP2(NN)=DXP(NN)*DXP(NN)+DYP(NN)*DYP(NN)+DZP(NN)*DZP(NN)-
     & DA(1)**2
      AD2M=MAX(AD2M,ABS(DP2(NN)))
      DXP(N1)=XX(NN1)-XX(NN3)
      DYP(N1)=YY(NN1)-YY(NN3)
      DZP(N1)=ZZ(NN1)-ZZ(NN3)
      DP2(N1)=DXP(N1)*DXP(N1)+DYP(N1)*DYP(N1)+DZP(N1)*DZP(N1)-
     & DA(2)**2
      AD2M=MAX(AD2M,ABS(DP2(N1)))
      DXP(N2)=XX(NN2)-XX(NN3)
      DYP(N2)=YY(NN2)-YY(NN3)
      DZP(N2)=ZZ(NN2)-ZZ(NN3)
      DP2(N2)=DXP(N2)*DXP(N2)+DYP(N2)*DYP(N2)+DZP(N2)*DZP(N2)-
     & DA(3)**2
      AD2M=MAX(AD2M,ABS(DP2(N2)))
c      write(*,*)' loop gamma',nn,ad2m
12    CONTINUE

c       write(*,*)'facto 2',facto1,facto2,facto3

c      write(*,*)'gamma1, esta iterando?',nn,iter,
c     &  dsqrt(dp2(nn)+da(1)**2),
c     &  dsqrt(dp2(n1)+da(2)**2),dsqrt(dp2(n2)+da(3)**2)

      IF(AD2M.GT.TOLL) THEN
      FLAG=.FALSE.
      
c      if(one.ne.two)goto 5666     
c      do i=1+natom,natom+nwat
c      i1=i+nwat
c      i2=i1+nwat
c      di1=(xx(i)-xx(i1))**2+(yy(i)-yy(i1))**2+(zz(i)-zz(i1))**2
c      di2=(xx(i)-xx(i2))**2+(yy(i)-yy(i2))**2+(zz(i)-zz(i2))**2
c      di3=(xx(i2)-xx(i1))**2+(yy(i2)-yy(i1))**2+(zz(i2)-zz(i1))**2
c      di1=dsqrt(di1)
c      di2=dsqrt(di2)
c      di3=dsqrt(di3)
c      write(*,*)'gamma2',i,i1,di1
c      write(*,*)'      ',i,i2,di2
c      write(*,*)'      ',i1,i2,di3
c      enddo

5666  continue

c       write(*,*)'facto 3',facto1,facto2,facto3

      NOFSET=NATOM
      DO 113 NN=1,NWAT
      N1=NN+NWAT
      N2=N1+NWAT
      NN1 = NOFSET + NN
      NN2 = NN1 + NWAT
      NN3 = NN2 + NWAT
      DP=DXP(NN)*DXT(NN)+DYP(NN)*DYT(NN)+DZP(NN)*DZT(NN)
      G=FACTO1*DP2(NN)/DP
c      write(*,*)' -1- ',dp,g,facto1,dp2(nn)

      XX(NN1) = XX(NN1) - G*DXT(NN) /WWMM(1)
      XX(NN2) = XX(NN2) + G*DXT(NN) /WWMM(2)
      YY(NN1) = YY(NN1) - G*DYT(NN) /WWMM(1)
      YY(NN2) = YY(NN2) + G*DYT(NN) /WWMM(2)
      ZZ(NN1) = ZZ(NN1) - G*DZT(NN) /WWMM(1)
      ZZ(NN2) = ZZ(NN2) + G*DZT(NN) /WWMM(2)

      DP=DXP(N1)*DXT(N1)+DYP(N1)*DYT(N1)+DZP(N1)*DZT(N1)
      G=FACTO2*DP2(N1)/DP
c      write(*,*)' -2- ',dp,g,facto2,dp2(n1)

      XX(NN1) = XX(NN1) - G*DXT(N1) /WWMM(1)
      XX(NN3) = XX(NN3) + G*DXT(N1) /WWMM(2)
      YY(NN1) = YY(NN1) - G*DYT(N1) /WWMM(1)
      YY(NN3) = YY(NN3) + G*DYT(N1) /WWMM(2)
      ZZ(NN1) = ZZ(NN1) - G*DZT(N1) /WWMM(1)
      ZZ(NN3) = ZZ(NN3) + G*DZT(N1) /WWMM(2)

      DP=DXP(N2)*DXT(N2)+DYP(N2)*DYT(N2)+DZP(N2)*DZT(N2)
      G=FACTO3*DP2(N2)/DP
c      write(*,*)' -3- ',dp,g,facto3,dp2(n2)

      XX(NN2) = XX(NN2) - G*DXT(N2) /WWMM(2)
      XX(NN3) = XX(NN3) + G*DXT(N2) /WWMM(2)
      YY(NN2) = YY(NN2) - G*DYT(N2) /WWMM(2)
      YY(NN3) = YY(NN3) + G*DYT(N2) /WWMM(2)
      ZZ(NN2) = ZZ(NN2) - G*DZT(N2) /WWMM(2)
      ZZ(NN3) = ZZ(NN3) + G*DZT(N2) /WWMM(2)
c      write(*,*)'durante',nn1,xx(nn1),yy(nn1),zz(nn1)
c      write(*,*)'durante',nn2,xx(nn2),yy(nn2),zz(nn2)
c      write(*,*)'durante',nn3,xx(nn3),yy(nn3),zz(nn3)
113   CONTINUE
c      write(*,*)'corrigio posic'

c       write(*,*)'facto 4',facto1,facto2,facto3
      ENDIF
c      do i=3,npart
c      write(*,*)' iter ',iter,xx(i),yy(i),zz(i),wwmm(1),wwmm(2),wwmm(3)
c      enddo

c       write(*,*)'facto 5',facto1,facto2,facto3
      IF(.NOT.FLAG) GOTO 611
      GO TO 1111

611   ITER=ITER+1 
c      write(*,*)'iter',iter
      IF(ITER.LT.MAXI5) GOTO 8000
      WRITE(6,1000) ITER,NN,AD2M,TOLL
c       write(6,*)abs(dp2(nn)),abs(dp2(n1)),abs(dp2(n2))
       write (6,*) 'problemas con los constraints de los solventes'
1000  FORMAT(' SLOW CONV,iter,NN,ad2m,  atol',2I10,2D20.10)
      IF (ITER.LT.MAXIT) GOTO 8000
      STOP
1111  CONTINUE

c       write(*,*)'facto 6',facto1,facto2,facto3
      DO I=1,NATOM
      GX(I)=ZERO
      GY(I)=ZERO
      GZ(I)=ZERO
      ENDDO
c      write(*,*)'puso a cero Gx(i) i=1,natom'

c       write(*,*)'facto 7',facto1,facto2,facto3
      NOFSET = NATOM
      DO 21 I=NATOM+1,NATSOL+NATOM
      IF(I.GT.NATOM+1)NOFSET=NOFSET+NNAT(I-1)
      DO 22 NN=NOFSET+1,NOFSET+NNAT(I)

      GX(NN)=XX(NN)-GX(NN)
      GY(NN)=YY(NN)-GY(NN)
      GZ(NN)=ZZ(NN)-GZ(NN)
c      write(*,*)'GX ',nn,gx(nn),gy(nn),gz(nn)
c      write(*,*)'gamma XX',nn,xx(nn),yy(nn),zz(nn)
    
c      FX(NN)=FX(NN)+GX(NN)*WWM(I)/DEL
c      FY(NN)=FY(NN)+GY(NN)*WWM(I)/DEL
c      FZ(NN)=FZ(NN)+GZ(NN)*WWM(I)/DEL

c      write(6,55)'gamm3 ',nn,gx(nn)*wwm(i)/del,gy(nn)*wwm(i)/del,
c     &   gz(nn)*wwm(i)/del
55    format(2x,'sale gam ',2i4,4x,3f16.8)

22    CONTINUE
21    CONTINUE

c      write(*,*)'Sale de gamma',itel

      RETURN
      END



