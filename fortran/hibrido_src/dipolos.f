      SUBROUTINE DIPOLOS(IT,NIN,NDIP,DV1,DV2,DV3,
     >   ux,uy,uz)
      INCLUDE 'param'
        INCLUDE 'COMM'
      DIMENSION DIP(150),DIPX(NAT),DIPY(NAT),DIPZ(NAT)
      DIMENSION XDIP(150),ZDIP(150)

      DBY1=2.54D0/0.52918D0/EE
      UDHH=ONE/DA(3)
      DOHH=DSQRT(DA(1)*DA(1)-(DA(3)*DA(3)/4.D0))
      DMH=DOHH-0.15D0
      UDMH=ONE/DMH

      DO I=1,NPART
      DIPX(I)=ZERO
      DIPY(I)=ZERO
      DIPZ(I)=ZERO
      ENDDO
       
      DIPT=ZERO
      DIPTX=ZERO
      DIPTY=ZERO
      DIPTZ=ZERO

      DO 5100 I=1+NATOM, NWAT+NATOM
      I1=I +NWAT
      I2=I1+NWAT

C-----X1 Y1 Z1: POSICION DEL SITIO M CON CARGA

      X1(I)=ALFA1*X(I)+ALFA2*(X(I1)+X(I2))
      Y1(I)=ALFA1*Y(I)+ALFA2*(Y(I1)+Y(I2))
      Z1(I)=ALFA1*Z(I)+ALFA2*(Z(I1)+Z(I2))

      DIPX(I)= (PC(I)*X1(I)+PC(I1)*X(I1)+PC(I2)*X(I2))*DBY1
      DIPY(I)= (PC(I)*Y1(I)+PC(I1)*Y(I1)+PC(I2)*Y(I2))*DBY1
      DIPZ(I)= (PC(I)*Z1(I)+PC(I1)*Z(I1)+PC(I2)*Z(I2))*DBY1

      DIP(I)=DSQRT( DIPX(I)*DIPX(I) + DIPY(I)*DIPY(I) +
     &              DIPZ(I)*DIPZ(I) )

      DIPTX = DIPTX + DIPX(I)
      DIPTY = DIPTY + DIPY(I)
      DIPTZ = DIPTZ + DIPZ(I)

C-----XXH YYH ZZH : COMP. VERSOR X(I) (HH)   
C-----XXX YYY ZZZZ: COMP. VERSOR Z(I) (OM)

      XXH=(X(I1)-X(I2))*UDHH
      YYH=(Y(I1)-Y(I2))*UDHH
      ZZH=(Z(I1)-Z(I2))*UDHH

      XXX=(0.5D0*(X(I1)+X(I2))-X1(I))*UDMH
      YYY=(0.5D0*(Y(I1)+Y(I2))-Y1(I))*UDMH
      ZZZZ=(0.5D0*(Z(I1)+Z(I2))-Z1(I))*UDMH

c---test1:
c      odX=dsqrt(xxh*xxh+yyh*yyh+zzh*zzh)
c      odZ=dsqrt(xxx*xxx+yyy*yyy+zzzz*zzzz)
c      write(*,*)'xxx yyy zzz',xxx,yyy,zzzz
c       w81=dsqrt((odx-1.D0)**2)
c       w82=dsqrt((odz-1.D0)**2)
c      if(w81.gt.(1.D-05))write(*,*)'odx.ne.one'
c      if(w82.gt.(1.D-05))write(*,*)'odz.ne.one'
c      write(*,*)'modX  modZ',i,odX,odZ

C-----MOMENTO DIPOLAR POR MOLECULA PROYECTADO EN X y Z
           
      XDIP(I)=DIPX(I)*XXH+DIPY(I)*YYH+DIPZ(I)*ZZH  
      ZDIP(I)=DIPX(I)*XXX+DIPY(I)*YYY+DIPZ(I)*ZZZZ  

c      WRITE(27,97)ITEL,I,XDIP(I),ZDIP(I)    

5100   CONTINUE    

C-----SUMO CONTRIBUCIONES DE LOS ATOMOS CUANTICOS(ux,uy,uz)
      DIPTX = DIPTX + UX
      DIPTY = DIPTY + UY
      DIPTZ = DIPTZ + UZ

C-----MOMENTO DIPOLAR TOTAL
      DIPT = DSQRT(DIPTX*DIPTX+DIPTY*DIPTY+DIPTZ*DIPTZ)

C-----PROYECCION SOBRE EJES DE INERCIA
      DV1= DIPTX*Xmi(1,1)+DIPTY*Xmi(2,1)+DIPTZ*Xmi(3,1)
      DV2= DIPTX*Xmi(1,2)+DIPTY*Xmi(2,2)+DIPTZ*Xmi(3,2)
      DV3= DIPTX*Xmi(1,3)+DIPTY*Xmi(2,3)+DIPTZ*Xmi(3,3)
   
      DIPTT=DSQRT(DV1*DV1+DV2*DV2+DV3*DV3)

      IF(MOD((IT-NIN),NDIP).EQ.0.AND.NDIP.GT.0)THEN
      WRITE(26,98)ITEL,(DIP(I),I=1+NATOM,NWAT+NATOM),DIPT
      WRITE(28,98)ITEL,DV1,DV2,DV3,DIPTT
      ENDIF


      RETURN
         
44    FORMAT(3X,I4,5x,4G14.7)      
97    FORMAT(I9,2X,I5,10G18.6)
98    FORMAT(I9,2X,10G18.6)
99    FORMAT(//,'  CARGAS  ',/)
91    FORMAT(2I5,3X,6G12.6)

      END

