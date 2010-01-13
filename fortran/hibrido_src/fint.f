      SUBROUTINE FINT(NT,IZ,NATSOL,FXH1,FYH1,FZH1)

C-----CALCULA ENERGIA Y FUERZAS COULOMBIANA
C-----Y LENN-JONES DE INTERACCION SIST. QUANTICO Y CLASICO
C-----(NUCLEOS)

      INCLUDE 'COMM'
      
      DIMENSION DX(NAT),DY(NAT),DZ(NAT),RIJSQ(NAT)
      DIMENSION DX1(NAT),DY1(NAT),DZ1(NAT),RIJSQ1(NAT)
      DIMENSION QIJ(NAT),IZ(NT),RIJ(NAT)
      
c     INTEGER US1,US2

      IF(NSPECQ.EQ.0) RETURN
      IF(NDFT.EQ.1)THEN   
      DO I=1,NATOM
      PC(I)=DBLE(IZ(I))*EE
C      write(45,*) 'carga',PC(I)/EE
      ENDDO
      ENDIF
     
C      DO I=1+NATOM,NPART
      DO I=1,NPART
      VCSS(I)=ZERO
      ENDDO
     

      DO 1554 I=1,NATOM
C---  Para que sea soluto soluto tambien.nano
c      do 133 j=I+1,NATOM
c      DX(j)= X(I) - X(j)
c
c      DY(j)= Y(I) - Y(j) 
c
c      DZ(j)= Z(I) - Z(j) 
c
c      DX1(j)=DX(j)
c      DY1(j)=DY(j)
c      DZ1(j)=DZ(j)
c      RIJSQ(j)= DX(j)*DX(j)+ DY1(j)* DY1(j)+DZ(j)*DZ(j)
c      QIJ(j)=PC(I)*PC(j)
c      RIJSQ1(j)=RIJSQ(j)
c133   continue

      DO 142 JJ=1,NWAT
      j=jj+natom

*-----DISTANCIAS ENTRE SITIOS CON CARGA: XYZ -----*


      IF (ICLSTR.EQ.1) THEN
      do ktr= 0,2
      koff=ktr*nwat
      DX(J+koff )=X(I)-X(J+koff)
      DY(J+koff)=Y(I)-Y(J+koff)
      DZ(J+koff)=Z(I)-Z(J+koff)


      RIJSQ(j+koff)=DX(J+koff)*
     > DX(J+koff)+DY(J+koff)
     > *DY(J+koff)+DZ(J+koff)
     > *DZ(J+koff)
      QIJ(J+koff)=PC(I)*PC(J+koff)
      enddo

      ELSE
*----------- Bulk ==> minima distancia -------------*
      do ktr= 0,2
      koff=ktr*nwat
 
      DX(J+koff)=(X(I)-X(J+koff))/BXLGTH
      DY(J+koff)=(Y(I)-Y(J+koff))/BXLGTH 
      DZ(J+koff)=(Z(I)-Z(J+koff))/BXLGTH 
 
      DX(J+koff)=(DX(J+koff)-ANINT(DX(J+koff)))*BXLGTH
      DY(J+koff)=(DY(J+koff)-ANINT(DY(J+koff)))*BXLGTH
      DZ(J+koff)=(DZ(J+koff)-ANINT(DZ(J+koff)))*BXLGTH


      RIJSQ(j+koff)=DX(J+koff)*
     > DX(J+koff)+DY(J+koff)
     > *DY(J+koff)+DZ(J+koff)
     > *DZ(J+koff)
      QIJ(J+koff)=PC(I)*PC(J+koff)
      enddo
      endif

*-----DISTANCIAS ENTRE SITIOS CON MASA: X1Y1Z1------*
      IF (ICLSTR.EQ.1) THEN
      do ktr= 0,2
      koff=ktr*nwat
      DX1(J+koff )=X1(I)-X1(J+koff)
      DY1(J+koff)=Y1(I)-Y1(J+koff)
      DZ1(J+koff)=Z1(I)-Z1(J+koff)
 
 
      RIJSQ1(j+koff)=DX1(J+koff)*
     > DX1(J+koff)+DY1(J+koff)
     > *DY1(J+koff)+DZ1(J+koff)
     > *DZ1(J+koff)
      enddo
 
      ELSE
*----------- Bulk ==> minima distancia -------------*
      do ktr= 0,2
      koff=ktr*nwat
 
      DX1(J+koff)=(X1(I)-X1(J+koff))/BXLGTH
      DY1(J+koff)=(Y1(I)-Y1(J+koff))/BXLGTH
      DZ1(J+koff)=(Z1(I)-Z1(J+koff))/BXLGTH
 
      DX1(J+koff)=(DX1(J+koff)-ANINT(DX1(J+koff)))*BXLGTH
      DY1(J+koff)=(DY1(J+koff)-ANINT(DY1(J+koff)))*BXLGTH
      DZ1(J+koff)=(DZ1(J+koff)-ANINT(DZ1(J+koff)))*BXLGTH
 
 
      RIJSQ1(j+koff)=DX1(J+koff)*
     > DX1(J+koff)+DY1(J+koff)
     > *DY1(J+koff)+DZ1(J+koff)
     > *DZ1(J+koff)
      enddo
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
         IF (RIJSQ1(J).LT.RCTSQ) THEN
        NC = NC +1
        JNFC(NC) = J
         ENDIF
1334     CONTINUE
 
        ENDIF


cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       Lennard Jones   soluto-solvente
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-------Entre Oxigenos clasicos y nucl. cuanticos  
C       Hay 1 solo sitio LJ en el solvente: 
C       oxigeno(nspecq+1)

        DO 1500 J = 1,NC

        SIG3  = PTFIVE * (SIGMA(I)+SIGMA(NSPECQ+1))
        EPSS  = DSQRT (EPS(I)*EPS(NSPECQ+1))
        IF (SIGMA(I).EQ.ZERO) THEN
        SIG3 = ZERO
        EPSS = ZERO
        ENDIF

c       write(*,*)'fint sig ',i,sig3,epss

        EPSI  = EPSS * FOUR * BOLTZF
        SIG3  = SIG3 * SIG3 * SIG3
        SIG6  = SIG3 * SIG3
        SIG12 = SIG6 * SIG6
        EE6   = EPSI * SIG6
        EE12  = EPSI * SIG12
        FF12  = 12.D0 * EE12
        FF6   = 6.D0  * EE6
        
c       write(*,*)'2fint sig ',i,sig6,sig12,ee6,ee12,ff6,ff12
       
        RRSQK1=ONE/RIJSQ1(JNFC(J))
        RRIJK1=DSQRT(RRSQK1)
        RR6=RRSQK1*RRSQK1*RRSQK1
        EGK=(EE12*RR6-EE6)*RR6
        FGK=(FF12*RR6-FF6)*RR6*RRSQK1
        VLJQC=VLJQC+EGK

888     format('1/r LJQC',2i4,10g14.6)


        AAX=DX1(JNFC(J))*FGK
        AAY=DY1(JNFC(J))*FGK
        AAZ=DZ1(JNFC(J))*FGK
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
       DO      JJ=0,2
               JJJ=JNFC(J) + JJ*NWAT
       
       RRSQK=ONE/RIJSQ(JJJ)
       RRIJK=DSQRT(RRSQK)        


       EGKC = QIJ(JJJ)*RRIJK
       FGKC = EGKC*RRSQK
       VCQC = VCQC + EGKC
       VCSS(JJJ)=VCSS(JJJ)+EGKC
       FQ(JJJ)=FQ(JJJ)-PC(I)*RRIJK
       
c        write(*,*)'COU', i,j,fgkc
c       write(*,889)i,j,1./rrijk,pc(i)/ee,pc(j)/ee,egkc
889    format('1/r VCQC ',2i4,10g14.6)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      AAX=DX(JJJ)*FGKC
      AAY=DY(JJJ)*FGKC
      AAZ=DZ(JJJ)*FGKC
C              write(6,*) 'AAX6',AAX        
 
C-----PARA LOS ATOMOS QUANTICOS (I) NO TRANSFIERE FZAS

      FX(I)=FX(I)+AAX
      FY(I)=FY(I)+AAY
      FZ(I)=FZ(I)+AAZ


C-----SI J ES UN OXIGENO (SOLVENTE)=> TRANSFIERE FZAS
      IF (JJ.EQ.0)THEN

      FX(JJJ)=FX(JJJ)-AAX*ALFA1
      FY(JJJ)=FY(JJJ)-AAY*ALFA1 
      FZ(JJJ)=FZ(JJJ)-AAZ*ALFA1 


      IIO=JJJ+NWAT
      FX(IIO)=FX(IIO)-AAX*ALFA2
      FY(IIO)=FY(IIO)-AAY*ALFA2
      FZ(IIO)=FZ(IIO)-AAZ*ALFA2

      IIO=IIO+NWAT
      FX(IIO)=FX(IIO)-AAX*ALFA2
      FY(IIO)=FY(IIO)-AAY*ALFA2
      FZ(IIO)=FZ(IIO)-AAZ*ALFA2


      ELSE

      FX(JJJ)=FX(JJJ)-AAX
      FY(JJJ)=FY(JJJ)-AAY 
      FZ(JJJ)=FZ(JJJ)-AAZ

      ENDIF
C      ENDIF
      enddo
1490  CONTINUE

 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

1554  CONTINUE


      if(UMBRE.EQ.1) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Muestreo paraguas I-K
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
 
       EGUM = (CKUM/2)*(RR-RUM)**2
       FGUM = -CKUM*(ONE - RUM/RR)
       VCUM = VCUM + EGUM
c       write(*,*) 'EGUM y FGUM',EGUM,FGUM
 
      AAX=DXU*FGUM
      AAY=DYU*FGUM
      AAZ=DZU*FGUM
 
      FX(US1)=FX(US1)+AAX
      FY(US1)=FY(US1)+AAY
      FZ(US1)=FZ(US1)+AAZ
 
      FX(US2)=FX(US2)-AAX
      FY(US2)=FY(US2)-AAY
      FZ(US2)=FZ(US2)-AAZ
 
C  144    continue
 
      elseif(UMBRE.EQ.2) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Muestreio Jarzynski
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DXU=X(US1)-X(US2)
      DYU=Y(US1)-Y(US2)
      DZU=Z(US1)-Z(US2)
c      write(*,*) J, DX(J), DY(J), DZ(J)


       RIJUMQ=DXU*DXU + DYU*DYU +DZU*DZU
       RR=DSQRT(RIJUMQ)
C       RRIJK=DSQRT(RRSQK)

       RUM2=RUM + ITEL*JVEL
       EGUM = (CKUM/2)*(RR-RUM2)**2
       FGUM = -CKUM*(ONE - RUM2/RR)
       VCUM = VCUM + EGUM
C       write(*,*) 'EGUM y FGUM',EGUM,FGUM,JVEL,RUM2

      AAX=DXU*FGUM
      AAY=DYU*FGUM
      AAZ=DZU*FGUM

      FX(US1)=FX(US1)+AAX
      FY(US1)=FY(US1)+AAY
      FZ(US1)=FZ(US1)+AAZ

      FX(US2)=FX(US2)-AAX
      FY(US2)=FY(US2)-AAY
      FZ(US2)=FZ(US2)-AAZ

       write(69,*) RUM2,FGUM 

c*****************dihedro!!!!

      elseif(UMBRE.EQ.3) THEN

       xmx=(y(1)-y(2))*(z(3)-z(2))-(z(1)-z(2))*(y(3)-y(2))

       xmy=-((x(1)-x(2))*(z(3)-z(2))-(z(1)-z(2))*(x(3)-x(2)))

       xmz=(x(1)-x(2))*(y(3)-y(2))-(y(1)-y(2))*(x(3)-x(2))

       xnx=(y(2)-y(3))*(z(4)-z(3))-(z(2)-z(3))*(y(4)-y(3))

       xny=-((x(2)-x(3))*(z(4)-z(3))-(z(2)-z(3))*(x(4)-x(3)))

       xnz=(x(2)-x(3))*(y(4)-y(3))-(y(2)-y(3))*(x(4)-x(3))

       sca=xmx*xnx+xmy*xny+xmz*xnz

       xm= xmx**2+xmy**2+xmz**2
       xm=dSQRT(xm)

       xn=xnx**2+xny**2+xnz**2
       xn=dSQRT(xn)

       if(xn*xm.eq.0.) pause 'dihedro problema puto!!!'
       arg=sca/(xm*xn)

       if(arg.ge.1.0) then 
       dih=0.
       pause
       elseif(arg.le.-1.0) then
       dih=PI
       pause
       else
       dih=dACOS(arg)
       endif
c       write(*,*) 'dihe',dih
       dih=dih*180/PI
c       write(*,*) 'dihe2',dih


        dtot = ckum*(dih-rum)
c        write(*,*) 'dtot',dtot,dih,rum,ckum
        prue=sca/(xn*xm)
        prue=(1.0-(prue)**2)
        if (prue.lt.1.0E-15.and.prue.gt.-1.0E-15) pause 'malo' 
        write(69,*) ITEL,dih

        prue=dsqrt(prue)
        dtot=dtot/prue

        do natdh=1,4
        do j=1,3
        idh=(natdh-1)*3+j

        dmx=0.0
        dmy=0.0
        dmz=0.0
        dnx=0.0
        dny=0.0
        dnz=0.0

        if(idh.eq.1) then
        dmy=Z(2)-Z(3)
        dmz=Y(3)-Y(2)
        elseif(idh.eq.4) then
        dmy=Z(3)-Z(1)
        dmz=Y(1)-Y(3)                
        dny=Z(3)-Z(4)                
        dnz=Y(4)-Y(3)                
        elseif(idh.eq.7) then
        dmy=Z(1)-Z(2)
        dmz=Y(2)-Y(1)
        dny=Z(4)-Z(2)
        dnz=Y(2)-Y(4)
        elseif(idh.eq.10) then
        dny=Z(2)-Z(3)
        dnz=Y(3)-Y(2)
        elseif(idh.eq.2) then
        dmx=Z(3)-Z(2)
        dmz=X(2)-X(3)
        elseif(idh.eq.5) then
        dmx=Z(1)-Z(3)
        dmz=X(3)-X(1)
        dnx=Z(4)-Z(3)
        dnz=X(3)-X(4)
        elseif(idh.eq.8) then
        dmx=Z(2)-Z(1)
        dmz=X(1)-X(2)
        dnx=Z(2)-Z(4)
        dnz=X(4)-X(2)
        elseif(idh.eq.11) then
        dnx=Z(3)-Z(2)
        dnz=X(2)-X(3)
        elseif(idh.eq.3) then
        dmx=Y(2)-Y(3)
        dmy=X(3)-X(2)
        elseif(idh.eq.6) then
        dmx=Y(3)-Y(1)
        dmy=X(1)-X(3)
        dnx=Y(3)-Y(4)
        dny=X(4)-X(3)
        elseif(idh.eq.9) then
        dmx=Y(1)-Y(2)
        dmy=X(2)-X(1)
        dnx=Y(4)-Y(2)
        dny=X(2)-X(4)
        elseif(idh.eq.12) then
        dnx=Y(2)-Y(3)
        dny=X(3)-X(2)
        endif


        dm=(xmx*dmx+xmy*dmy+xmz*dmz)/xm
        dn=(xnx*dnx+xny*dny+xnz*dnz)/xn
        dmn=xm*dn+xn*dm
       dscalar=xnx*dmx+xmx*dnx+xny*dmy+xmy*dny+xnz*dmz+xmz*dnz

c         write(66,*)'dmn dn dm dmx dmy dmz dny dnx dnz'
c         write(66,*) dmn, dn, dm, dmx, dmy, dmz, dny, dnx, dnz
        if(j.eq.1) then
         FX(natdh)=FX(natdh)+dtot*(dscalar*xm*xn-dmn*sca)
     >  /(xn*xm)**2
        elseif(j.eq.2) then
         FY(natdh)=FY(natdh)+dtot*(dscalar*xm*xn-dmn*sca)
     >  /(xn*xm)**2

c        write(66,*) natdh,dtot,(dscalar*xm*xn-dmn*sca),(xn*xm)**2

        elseif(j.eq.3) then
         FZ(natdh)=FZ(natdh)+dtot*(dscalar*xm*xn-dmn*sca)
     >  /(xn*xm)**2
         endif
        enddo
        enddo

          VCUM = VCUM + (ckum*(dhi-rum)**2)/2

      ENDIF




C    DX1(J)=X1(I)-X1(J)
C     DY1(J)=Y1(I)-Y1(J)
C     DZ1(J)=Z1(I)-Z1(J)
c     RIJ(J)=DX1(J)*DX1(J)+DY1(J)*DY1(J)+DZ1(J)*DZ1(J)
c     RIJ(J)=DSQRT(RIJ(J))


C-----G(R) ENTRE ATOMOS CUANTICOS
c      IRR = NINT(RIJ(J)/DELR)
C      IRR = MIN(IGGRID,IRR)
C      KG(M1,IRR)=KG(M1,IRR)+1
      
C80    CONTINUE

cC--------------------------------------
C-----PARA EL DIMERO DE AGUA :G(r) O-O
c      do i=1,npart
c      write(6,233)i,x(i),y(i),z(i)
c      enddo
233   format(i5,3f10.5)

c      DDX=X(1)-X(4)
c      DDY=Y(1)-Y(4)
c      DDZ=Z(1)-Z(4)
c      RRIJ=DSQRT(DDX*DDX+DDY*DDY+DDZ*DDZ)
c       write(*,*)'fint : ',rrij
c      IRR = NINT(RRIJ/DELR)
c      IRR = MIN(IGGRID,IRR)
c      KG(4,IRR)=KG(4,IRR)+1
C      IF(NDFT.EQ.1)THEN
c      DO I=1,NATOM
cC      PC(I)=IZ(I)
cc      ENDDO
c      ENDIF

ccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
cc-----Control de fzas en fint
      fxx=zero
      fyy=zero
      fzz=zero
      do i=1,npart
      fxx=fxx+fx(i)
      fyy=fyy+fy(i)
      fzz=fzz+fz(i)
      enddo
     

c      write(*,*)'hx hy hz ',dsqrt(hxx**2+hyy**2+hzz**2)
      if(dsqrt(fxx**2+fyy**2+fzz**2).gt.1.d-04)then
      write(*,*)'FZA TOTAL NE ZERO EN FINT  '  
      write(*,78)itel,fxx,fyy,fzz
      endif
78    format(i9,3g15.7)

      RETURN
      END

      

