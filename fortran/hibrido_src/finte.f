      SUBROUTINE FINTE(IZ,NATSOL,FXH1,FYH1,FZH1)


C-----CALCULA ENERGIA Y FUERZAS COULOMBIANA
C-----Y LENN-JONES DE INTERACCION SIST. QUANTICO Y CARGAS 
C-----CLASICAS DEL CAPACITOR. LAS CARGAS DEL CAPACITOR SON
C-----ATOMOS DE OXIGENO CON CARGA DEPENDIENTE DEL CAMPO

      INCLUDE 'param'
        INCLUDE 'COMM'

      DIMENSION DX(NAT),DY(NAT),DZ(NAT),RIJSQ(NAT)
      DIMENSION QIJ(NAT),IZ(NT),RIJ(NAT)

      IF(NSPECQ.EQ.0) RETURN
      IF(NDFT.EQ.1)THEN
      DO I=1,NATOM
      PC(I)=DBLE(IZ(I))*EE
      ENDDO
      ENDIF

      DO I=1,NPART
      VCSS(I)=ZERO
      ENDDO
      
      DO 1554 I=1,NATOM
      DO 142 J=NATOM+1,NPART

*-----DISTANCIAS ENTRE SITIOS CON MASA: XYZ------*
      IF (ICLSTR.EQ.1) THEN
      DX(J)=X(I)-X(J)
      DY(J)=Y(I)-Y(J)
      DZ(J)=Z(I)-Z(J)

      RIJSQ(j)=DX(J)*
     > DX(J)+DY(J)
     > *DY(J)+DZ(J)
     > *DZ(J)
      QIJ(J)=PC(I)*PC(J)
      ELSE
*----------- Bulk ==> minima distancia -------------*

      DX(J)=(X(I)-X(J))/BXLGTH
      DY(J)=(Y(I)-Y(J))/BXLGTH
      DZ(J)=(Z(I)-Z(J))/BXLGTH

      DX(J)=(DX(J)-ANINT(DX(J)))*BXLGTH
      DY(J)=(DY(J)-ANINT(DY(J)))*BXLGTH
      DZ(J)=(DZ(J)-ANINT(DZ(J)))*BXLGTH


      RIJSQ(j)=DX(J)*
     > DX(J)+DY(J)
     > *DY(J)+DZ(J)
     > *DZ(J)
      endif

142   CONTINUE

        IF (ICLSTR.EQ.1) THEN
        NC = 0
        DO 9133 J = natom+1,NWAT+NATOM
        NC = NC +1
        JNFC(NC)=J
9133    CONTINUE

        ELSE

        NC = 0
        DO 1334 J = natom+1,NWAT+NATOM
         IF (RIJSQ(J).LT.RCTSQ) THEN
        NC = NC +1
        JNFC(NC) = J
         ENDIF
1334     CONTINUE

        ENDIF

cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       Lennard Jones   soluto-solvente
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-------Entre Oxigenos clasicos y nucl. cuanticos
C       Hay 1 solo sitio LJ 
C       oxigeno(nspecq+1)

        DO 1500 J = 1,NC

        SIG3  = PTFIVE * (SIGMA(I)+SIGMA(NSPECQ+1))
        EPSS  = DSQRT (EPS(I)*EPS(NSPECQ+1))
        IF (SIGMA(I).EQ.ZERO) THEN
        SIG3 = ZERO
        EPSS = ZERO
        ENDIF

        EPSI  = EPSS * FOUR * BOLTZF
        SIG3  = SIG3 * SIG3 * SIG3
        SIG6  = SIG3 * SIG3
        SIG12 = SIG6 * SIG6
        EE6   = EPSI * SIG6
        EE12  = EPSI * SIG12
        FF12  = 12.D0 * EE12
        FF6   = 6.D0  * EE6

        RRSQK=ONE/RIJSQ(JNFC(J))
        RRIJK=DSQRT(RRSQK)
c       WRITE(*,*) RRIJK
        RR6=RRSQK*RRSQK*RRSQK
        EGK=(EE12*RR6-EE6)*RR6
        FGK=(FF12*RR6-FF6)*RR6*RRSQK
        VLJQC=VLJQC+EGK

888     format('1/r LJQC',2i4,10g14.6)

        AAX=DX(JNFC(J))*FGK
        AAY=DY(JNFC(J))*FGK
        AAZ=DZ(JNFC(J))*FGK
c       write(6,*) 'AAX2',AAX

        FX(I)=FX(I)+AAX
        FY(I)=FY(I)+AAY
        FZ(I)=FZ(I)+AAZ

        FX(JNFC(J))=FX(JNFC(J))-AAX
        FY(JNFC(J))=FY(JNFC(J))-AAY
        FZ(JNFC(J))=FZ(JNFC(J))-AAZ

C        ENDIF
1500    CONTINUE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Coulomb soluto-solvente
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       DO 1490 J=1,NC
               JJJ=JNFC(J) 

       RRSQK=ONE/RIJSQ(JJJ)
       RRIJK=DSQRT(RRSQK)
       EGKC = QIJ(JJJ)*RRIJK
c      WRITE(*,*) 'QIJ(JJJ)',QIJ(JJJ)
       FGKC = EGKC*RRSQK
       VCQC = VCQC + EGKC
       VCSS(JJJ)=VCSS(JJJ)+EGKC
       FQ(JJJ)=FQ(JJJ)-PC(I)*RRIJK

889    format('1/r VCQC ',2i4,10g14.6)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      AAX=DX(JJJ)*FGKC
      AAY=DY(JJJ)*FGKC
      AAZ=DZ(JJJ)*FGKC


      FX(I)=FX(I)+AAX
      FY(I)=FY(I)+AAY
      FZ(I)=FZ(I)+AAZ


      FX(JJJ)=FX(JJJ)-AAX
      FY(JJJ)=FY(JJJ)-AAY
      FZ(JJJ)=FZ(JJJ)-AAZ


1490  CONTINUE

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

1554  CONTINUE

      if(UMBRE.EQ.1) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Muestreo d1-d2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DXU=X(US1)-X(US2)
      DYU=Y(US1)-Y(US2)
      DZU=Z(US1)-Z(US2)
c      write(*,*) J, DX(J), DY(J), DZ(J)


       RIJUMQ=DXU*DXU + DYU*DYU +DZU*DZU
       RR=DSQRT(RIJUMQ)
C       RRIJK=DSQRT(RRSQK)
c       write(*,*) 'RIJUMQ',RIJUMQ,RR,RUM

       write(69,*) ITEL,RR
       CALL FLUSH(69)
       endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-------------------------------------------------------------------------
      fxx=zero
      fyy=zero
      fzz=zero
      do i=1,npart
      fxx=fxx+fx(i)
      fyy=fyy+fy(i)
      fzz=fzz+fz(i)
      enddo


      if(dsqrt(fxx**2+fyy**2+fzz**2).gt.1.d-04)then
      write(*,*)'FZA TOTAL NE ZERO EN FINT  '
      write(*,78)itel,fxx,fyy,fzz
      endif
78    format(i9,3g15.7)

      RETURN
      END

