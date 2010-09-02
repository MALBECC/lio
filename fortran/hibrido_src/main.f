C
C   PROGRAMA PARA MOVER UN SISTEMA FORMADO POR 2 SUBSISTEMAS:
C   SOLVENTE==> TIP4P o TIP4P-FQ-WATER CON Q-FLUCTUANTES
C   SOLUTO  ==> SISTEMA QUANTICO (DFT)
C   
C
C
C-----DIMENSIONA PARAMETROS SIST.QUANTICO:
C-----(init.f de Dario)
c
c MAIN SUBROUTINE --------------------------------------
C DFT calculation with gaussian basis sets
c-------------------------------------------------------
c
c
c PARAMETERS - DYNAMICAL VECTOR ONLY -------------------
c
c ngDyn : number of atoms * number basis functions
c ngdDyn: number of atoms * number auxiliar functions
c Ngrid : number of grid points (LS-SCF part)
c norbit : number of MO
c
c Ngrid may be set to 0 , in the case of using Num. Integ.

      INCLUDE 'COMM'
      INCLUDE 'param'
      INTEGER SPC
      COMMON /tipsol/SPC
      parameter (ngDyn=700)
      parameter (ngdDyn=700)
c      parameter (ngDyn=450)
c      parameter (ngdDyn=450)
c      parameter (norbit=250,Ngrid=8000)
      parameter (norbit=80,Ngrid=6000)

      parameter (ng3=4*ngDyn)
c     parameter (ng2=7*ngDyn*(ngDyn+1)/2+3*ngdDyn*(ngdDyn+1)/2+
c    >           ngDyn+ngDyn*norbit+Ngrid)
c---- para version en memoria
      parameter (ng2=5*ngDyn*(ngDyn+1)/2+3*ngdDyn*(ngdDyn+1)/2+
     >           ngDyn+ngDyn*norbit+Ngrid+ngDyn*(NgDyn+1)/2*NgdDyn)

      dimension XW(ngDyn,ng3),XXW(ngdDyn,ngdDyn),RMM(ng2)

C-----DIMENSIONES DE TODO
      DIMENSION AV(500)
      logical NORM,ATRHO,VCINP,DIRECT,EXTR,dens,write1
      logical OPEN,SVD,SHFT,GRAD,BSSE,integ,field,sol,free
      logical exter,MEMO
      dimension r(nt,3),nshelld(0:3),nshell(0:3)
      dimension cd(ngd,nl),ad(ngd,nl),Nucd(ngd),ncontd(ngd)
      dimension c(ng,nl),a(ng,nl),Nuc(ng),ncont(ng),Iz(nt)
      dimension Em(ntq+nss),Rm(ntq+nss),alpha(nss),pcx(nt)
      dimension xte(nat),yte(nat),zte(nat)
      dimension fxte(nat),fyte(nat),fzte(nat) 
      dimension f1(nt,3),Pm(nt),q(ntq)
      DIMENSION IZTE(NAT)
      DIMENSION HISTO(100,100),HISTO2(100)
      CHARACTER*4 date

      common /sol1/Nsol,natsol,alpha,Em,Rm,sol,free,pcx

      common /dyn/Pm
      common /fit/ Nang_,dens,integ,Iexch,igrid,igrid2

#ifdef G2G
      call g2g_init()
#endif

      TEMPAV=ZERO
      TEMPOL=ZERO
      TEMPSLV=ZERO
      TEMPSLT=ZERO
      IFORT=70


C--------------------------------------------------------
C--------------------------------------------------------

C-----LLAMA A 'INICIO':LEE TODO SOBRE EL SISTEMA CLASICO
      CALL INICIO(NATSOL,NDIP,IDIPCOR,PMAX,PZMAX)

      KBER = BER1
      IPINPUT0 =  ITERM

C-----LLAMA A 'DRIVE' :LEE TODO SOBRE EL SISTEMA QUANTICO
C-----(AUNQUE ICON=1 O -1 LEE DE DRIVE, SI ICON.NE.0 LEE LAS 
C-----POSICIONES PERO MAS ADELANTE LLAMA A CONFIG Y REESCRIBE 
C-----LAS POSICIONES NUCLEARES.
c      IF(NSPECQ.NE.0)THEN
      CALL DRIVE(ng2,ngDyn,ngdDyn,RMM,XW,XXW,MEMO,NORM,natom,Iz,r,Nuc,
     & M,ncont,nshell,c,a,Nucd,Md,ncontd,nshelld,cd,ad,E,
     & nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write1)
c      ENDIF

       date='date'
       CALL SYSTEM(date)
       write(*,*)
     
C-----LLAMA A 'CORECT': CALCULA TODO LO QUE NECESITA DESPUES

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(ELFIELD.EQ.1) THEN
       CALL CORECTE(NATSOL,RM,EM,NTQ,NSS,PM,NT,AVNPU,
     &  PMAX,PZMAX,HISTO,HISTO2,DELP,DELPZ)
       GOTO 6778
       ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   
      IF(SPC.EQ.1)THEN 
      CALL CORECT2(NATSOL,RM,EM,NTQ,NSS,PM,NT,AVNPU,
     &  PMAX,PZMAX,HISTO,HISTO2,DELP,DELPZ)
      ELSE
      CALL CORECT(NATSOL,RM,EM,NTQ,NSS,PM,NT,AVNPU,
     &  PMAX,PZMAX,HISTO,HISTO2,DELP,DELPZ)
      ENDIF 
 6778 CONTINUE
     
C-----TRANSFORMA NUMERACION Y UNIDADES DARIO A 
C-----LAS MIAS (COORDS)
      IF(ABS(ICON).EQ.0)THEN

      DO I=1+NATOM,NPART
      IZTE(I)=IZ(I)
      ENDDO

      DO I=1,NATOM
      X(I)=R(I,1)*A0
      Y(I)=R(I,2)*A0
      Z(I)=R(I,3)*A0
      ENDDO

C===========================================
      IF(ELFIELD.EQ.1) THEN
      
      DO I=NATOM+1,NATOM+NWAT
      X(I)=R(I,1)*A0
      Y(I)=R(I,2)*A0
      Z(I)=R(I,3)*A0
      IZ(I)=IZTE(I)
      ENDDO

      GOTO 6779
      ENDIF
      
C===========================================      

      II=NATOM
      DO I=NATOM+1,NATOM+NWAT*NATSOL,3
      II=II+1
      X(II)=R(I,1)*A0
      Y(II)=R(I,2)*A0
      Z(II)=R(I,3)*A0
      IZ(II)=IZTE(I)
      ENDDO
         
      DO I=NATOM+1,NATOM+NWAT*NATSOL,3
      II=II+1
      X(II)=R(I+1,1)*A0
      Y(II)=R(I+1,2)*A0
      Z(II)=R(I+1,3)*A0
      IZ(II)=IZTE(I+1)
      ENDDO

      DO I=NATOM+1,NATOM+NWAT*NATSOL,3
      II=II+1
      X(II)=R(I+2,1)*A0
      Y(II)=R(I+2,2)*A0
      Z(II)=R(I+2,3)*A0
      IZ(II)=IZTE(I+2)
      ENDDO

 6779 CONTINUE

C-----SI COMIENZA DEL ORDEN (ICON=0) LEYO CONFIGURACION DE
C-----DARIO Y EMPIEZA CON VELOCIDAD NULA
C-----(EMPEZAR DEL ORDEN QUIERE DECIR DE DRIVE)
      DO I=1,NPART
      X0(I)=X(I)
      Y0(I)=Y(I)
      Z0(I)=Z(I)
      VX(I)=ZERO
      VY(I)=ZERO
      VZ(I)=ZERO
      VX0(I)=ZERO
      VY0(I)=ZERO
      VZ0(I)=ZERO
      VQ(I)=ZERO
      VQ0(I)=ZERO
      ENDDO

      ENDIF

C-----CONTROL: COHERENCIA ENTRE ENTRADA DARIO Y MIA
c      IF(NSPECQ.EQ.0)NSOL=NWAT
      IF(NWAT.NE.NSOL)THEN
      WRITE(*,*)'NSOL DEBE SER = A NWAT',NSOL,NWAT
      PAUSE
      ENDIF
      IF(NPART.NE.(NATOM+NSOL*NATSOL))THEN
      WRITE(*,*)'ALGO INCOHERENTE NPART,NT',NPART,(NATOM+NSOL*NATSOL)
      STOP
      ENDIF
     
C-----PONE A CERO LAS CANTIDADES CANONICAS
      S0 = ZERO
      ST = ZERO
      SD = ZERO
      SS = ZERO
      SD0= ZERO
     
      S0Q = ZERO
      STQ = ZERO
      SDQ = ZERO
      SSQ = ZERO
      SD0Q= ZERO

C-----LLAMA A 'CONFIG': COMO EMPEZAR (NUCLEOS)
      CALL CONFIG (ANG,NATSOL,ntq,q,HISTO,HISTO2)


C-----ARREGLA LA NUMERACION PARA EL IZ(I)
C-----IZ(I):NUMERO ATOMICO
      DO I=1,NPART
      IZTE(I)=IZ(I)
      ENDDO

      DO II=1,NATSOL
      II1=II-1
      DO I=NATOM+1,NATOM+NWAT
      K=(I-(NATOM+1))*(NATSOL-1)  
      IZ(I+II1*NWAT)=IZTE(I+II1+K)
      ENDDO
      ENDDO


C-----SI QUIERO PONER LAS QS INICIALES 
C-----DEBE SER IZZ=1

      IF (IZZ.EQ.1)THEN

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(ELFIELD.EQ.1) THEN

      DO I=1+NATOM,0.5*NWAT+NATOM
      II= I+0.5*NWAT

      PC(I) = ZZZ(1)*EE
      IZ(I) = 8
      PC(II)= -ZZZ(1)*EE
      IZ(II) = 7
      ENDDO
      
c     DO I=1+NATOM,NWAT+NATOM
c     WRITE(*,*) I,PC(I)
c     ENDDO
          
      GOTO 1090
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DO I=1+NATOM,NWAT+NATOM
      I1= I+NWAT
      I2= I+2*NWAT
      
      PC(I) = ZZZ(1)*EE
      PC(I1)= ZZZ(2)*EE
      PC(I2)= ZZZ(3)*EE
      PC0(I)  = PC(I)
      PC0(I1) = PC(I1)
      PC0(I2) = PC(I2)
 
      ENDDO
 1090 CONTINUE 
      ENDIF


C-----SI ICON=0 PERO QUIERO LEER LAS CARGAS DEL SOLVENTE(IZZ=2)
C-----(LAS LEE CON MI NUMERACION Y UNID.DE E)
      IF(ICON.EQ.0.AND.IZZ.EQ.2)THEN
      OPEN(74,FILE='cargas')
      READ(74,*) NATOM,NPART
      DO I=NATOM+1,NPART
      READ(74,*)NU,PC(NU)
      PC(NU)=PC(NU)*EE
      PC0(NU)=PC(NU)
      ENDDO
      ENDIF


      REWIND 8
      NIN = ITEL +2
      NOUT = NIN + NTIME
      GGRID = 1000.D0
      IGGRID=INT(GGRID)
      GRIL = 100.D0
      IGRIL= INT(GRIL)
      DELR2= BXLGTH/GRIL
      DELR = BXLGTH/GGRID
 

C-----GENERA UNA FOTO (ini.xyz)INICIAL
        WRITE(12,50) NPART
        WRITE(12,*)
        DO I=1,NPART
        IF(NDFT.EQ.1)THEN     

         IF(I.LE.NATOM)THEN
          WRITE(12,101)   AT(I),X(I),Y(I),Z(I)
         ELSE
          WRITE(12,101)   AT(I),X(I),Y(I),Z(I),PC(I)/EE
         ENDIF

        ELSE
         WRITE(12,101)   AT(I),X(I),Y(I),Z(I),PC(I)/EE
        ENDIF

        ENDDO
        WRITE(12,*)'inicio '
        

C-- Para reacomodar el solvente a TIP4P o SPC:  
       do i=1,npart
       xx(i)=x(i)
       yy(i)=y(i)
       zz(i)=z(i)
       enddo
       CALL GAMMA(NATSOL)
       do i=1,npart
       x(i)=xx(i)
       y(i)=yy(i)
       z(i)=zz(i)
       enddo
C-- Fin reacomodo 


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
CC   AQUI COMIENZA EL LOOP
CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DO 390 IT = NIN , NOUT
      ITEL = IT -1


C-----LLAMA A 'NEWXYZ':PONE EN X1 Y1 Z1 SITIOS REALES
C-----Y EN XYZ LOS SITIOS CON CARGA (TIP4P)  
c     write(*,*) 'ANTES DE NEW', X(6),Y(6),Z(6)   
      CALL NEWXYZ(NATSOL,F1,NT)
c     write(*,*) 'DESPUES DE NEW',X(6),Y(6),Z(6)
   
      Vcerca = ZERO
      Vlejos = ZERO 
      VCT    = ZERO
      VLJ    = ZERO
      VOO    = ZERO
      VCOO   = ZERO
      VOH    = ZERO
      VHH    = ZERO
      VTOT   = ZERO
      VSELF  = ZERO
      VLJQC  = ZERO
      VCQC   = ZERO
      VCSSO  = ZERO
      VCSSH  = ZERO
      VCSST  = ZERO
      FXH1   = ZERO
      FYH1   = ZERO
      FZH1   = ZERO
      FXH11  = ZERO
      FYH11  = ZERO
      FZH11  = ZERO
      FXH2   = ZERO
      FYH2   = ZERO
      FZH2   = ZERO
      
       
 
 
      IF(NWAT.EQ.0)GOTO 575

C-----LLAMA A 'FSPHER':
c     write(*,*)'SPC:',SPC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(ELFIELD.EQ.1) THEN
      CALL FSPHER(NATSOL)
      GOTO 6880
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF(SPC.EQ.1) THEN
      CALL FSPHER2(NATSOL)
c     write(*,*)'fspher2'
      ELSE
      CALL FSPHER(NATSOL)
c     write(*,*)'fspher'
      ENDIF
 
 6880 CONTINUE

575   CONTINUE

C-----LLAMA A 'FINT': INT. SOLV-SLT(NUCLEOS)
C      write(*,*)'a fint',FX(1),FY(1),Fz(1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(ELFIELD.EQ.1) THEN
      CALL FINTE(NT,IZ,NATSOL,FXH1,FYH1,FZH1)
      
c     WRITE(*,*) VLJQC,VCQC
      GOTO 6881
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF(SPC.EQ.1) THEN
      CALL FINT2(NT,IZ,NATSOL,FXH1,FYH1,FZH1)
      ELSE
      CALL FINT(NT,IZ,NATSOL,FXH1,FYH1,FZH1)
      ENDIF 
C      write(*,*)' fint',FX(1),Fy(1),Fz(1)

 6881 CONTINUE

C575   CONTINUE
C-----VUELVE A PONER EN XYZ LOS SITIOS REALES
C-----Y EN X1 Y1 Z1 LOS SITIOS CON CARGA (TIP4P)
       DO I = 1, NPART

       XT = X1(I)
       YT = Y1(I)
       ZT = Z1(I)

       X1(I) = X(I)
       Y1(I) = Y(I)
       Z1(I) = Z(I)

       X(I) = XT
       Y(I) = YT
       Z(I) = ZT
       ENDDO

        IF(NDFT.NE.1)GOTO 9000
      IF (NPAS.GT.0.AND.MOD(ITEL,NPAS).EQ.0) THEN

C-----LLAMA A 'DODA':PASA DE MI NUMERACION A LA DE DARIO (COORDS)

      CALL DODA(NATSOL,R,NT,IZ)

C-----LLAMA A 'SCF': CALCULA ENERGIA SUBS. QUANTICO
      IF(NPAS.NE.1)THEN
       WRITE(*,*)
       WRITE(*,*)'DFT-   STEP No: ',ITEL
       WRITE(*,*)
      ENDIF


#ifdef G2G
      write(*,*) 'carga de posiciones, igrid2'
			call g2g_reload_atom_positions(igrid2)
#endif
         IF(OPEN)THEN

      CALL SCFop(MEMO,NORM,NATOM,IZ,R,NUC,M,NCONT,NSHELL,C,A
     >    ,NUCD,MD,NCONTD,NSHELLD,CD,AD
     >    ,RMM,XW,XXW,E,
     >    NOPT,PC,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,NUNP,GOLD,
     >    TOLD,WRITE1,FQQ,Q,
     >    IT,ITEL,NIN,IPR1,E1s,EAC,
     >    ux,uy,uz,NPAS)

          ELSE

      CALL SCF(MEMO,NORM,NATOM,IZ,R,NUC,M,NCONT,NSHELL,C,A
     >    ,NUCD,MD,NCONTD,NSHELLD,CD,AD
     >    ,RMM,XW,XXW,E,
     >    NOPT,PC,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,NUNP,GOLD,
     >    TOLD,WRITE1,FQQ,Q,
     >    IT,ITEL,NIN,IPR1,E1s,EAC,
     >    ux,uy,uz,NPAS)

          ENDIF

      ESCF = E
      E1s0 = E1s
      EKS  = ESCF-E1s0

#ifdef G2G
			write(*,*) 'cambio de grilla para fuerza (igrid)'
      if (igrid.ne.igrid2) then
        call g2g_new_grid(igrid)
      endif
#endif

C-----LLAMA A RUTINAS QUE CALCULAN GRADIENTES DE EKS:

       
C---  INT1G:GRADIENTES DE INTEGRALES DE 1E (Q)
      CALL INT1G(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
     > c,a,RMM,En,f1,FXH2,FYH2,FZH2)
      

C---  INT3G:GRADIENTES DE INTEGRALES DE 2E (Q)
      CALL INT3G(NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     > Nucd,Md,ncontd,nshelld,cd,ad,RMM,E2,f1,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,
     > told,write,FXH2,FYH2,FZH2)


C---  INTSG:GRADIENTES DE INTEGRAL DE OVERLAP (Q)
      CALL INTSG(NORM,natom,nsol,natsol,r,Nuc,M,Md,ncont,
     > nshell,c,a,RMM,f1,FXH2,FYH2,FZH2)
       FXH2 =  f1(1,1)
       FYH2 =  f1(2,1)
       FZH2 =  f1(3,1)

543   CONTINUE     

C-----INTSOLG:GRADIENTES DE INTEGRAL RHO*V (Q-C)
      CALL INTSOLG(NORM,natom,Nsol,natsol,r,Nuc,Iz,M,Md,
     > ncont,nshell,c,a,pc,RMM,f1,FXH11,FYH11,FZH11)
       
C-----UNIDADES DE DARIO POR LAS MIAS PARA LAS FZAS SOBRE SLT
      DO I=1,NATOM
      FX(I)=FX(I)-F1(I,1)/A0*HH
      FY(I)=FY(I)-F1(I,2)/A0*HH
      FZ(I)=FZ(I)-F1(I,3)/A0*HH
      ENDDO
c      write(*,*) 'd intsolG',FX(1),FY(1),FZ(1)
C-----ACOPLAMIENTO COULOMBIANO SOLUTO-SOLVENTE POR ATOMO
C-----A PARTIR DE ACA, VCSS ES LA SUMA DE LAS INTERACCIONES
C-----1) ELECTRONES CON NUCLEOS CLASICOS (EAC(I))
C-----2) NUCLEOS CLASICOS CON NUCLEOS CUANTICOS (VCSS(I))

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(ELFIELD.EQ.1) THEN
      I1=NATOM
      DO I=1+NATOM,NPART
      I1=I1+1
      VCSS(I1)=VCSS(I1)*FACTOR + EAC(I)*HH*FACTOR
      VCSSO = VCSSO + VCSS(I1)
      ENDDO
      GOTO 6882
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      I1=NATOM  ! O
      DO I=1+NATOM,NPART,3
      I1=I1+1
      VCSS(I1)=VCSS(I1)*FACTOR + EAC(I)*HH*FACTOR
      VCSSO = VCSSO + VCSS(I1)
      ENDDO
 
      DO I=1+NATOM,NPART,3
      I1=I1+1   ! H
      VCSS(I1)=VCSS(I1)*FACTOR + EAC(I+1)*HH*FACTOR
      VCSSH = VCSSH + VCSS(I1)
      ENDDO
 
      DO I=1+NATOM,NPART,3
      I1=I1+1   ! H
      VCSS(I1)=VCSS(I1)*FACTOR + EAC(I+2)*HH*FACTOR
      VCSSH = VCSSH + VCSS(I1)
      ENDDO

 6882 CONTINUE

      VCSST = VCSST + VCSSO + VCSSH
c     WRITE(*,*) 'VCSST VCSSO VCSSH',VCSST,VCSSO,VCSSH
C-----UNIDADES DE DARIO A LAS MIAS PARA LAS FZAS SOBRE SOLV.

C=============== ELFIELD=1 ==========================
 
      IF(ELFIELD.EQ.1) THEN
      NN=NATOM
      DO N=NATOM+1,NATOM+NWAT
      NN=NN+1
      FX(NN)=FX(NN)-F1(N,1)/A0*HH
      FY(NN)=FY(NN)-F1(N,2)/A0*HH
      FZ(NN)=FZ(NN)-F1(N,3)/A0*HH
      FQ(NN)=FQ(NN)+FQQ(N)/EE*HH
      enddo
      GOTO 6883
      ENDIF
C====================================================

      NN=NATOM
      DO N=NATOM+1,NATOM+NWAT*3,3
      NN=NN+1
      FX(NN)=FX(NN)-F1(N,1)/A0*HH
      FY(NN)=FY(NN)-F1(N,2)/A0*HH
      FZ(NN)=FZ(NN)-F1(N,3)/A0*HH
      FQ(NN)=FQ(NN)+FQQ(N)/EE*HH
      enddo

      DO N=NATOM+1,NATOM+NWAT*NATSOL,3
      NN=NN+1
      FX(NN)=FX(NN)-F1(N+1,1)/A0*HH
      FY(NN)=FY(NN)-F1(N+1,2)/A0*HH
      FZ(NN)=FZ(NN)-F1(N+1,3)/A0*HH
      FQ(NN)=FQ(NN)+FQQ(N+1)/EE*HH
      ENDDO
 
      DO N=NATOM+1,NATOM+NWAT*NATSOL,3
      NN=NN+1
      FX(NN)=FX(NN)-F1(N+2,1)/A0*HH
      FY(NN)=FY(NN)-F1(N+2,2)/A0*HH
      FZ(NN)=FZ(NN)-F1(N+2,3)/A0*HH
      FQ(NN)=FQ(NN)+FQQ(N+2)/EE*HH
      ENDDO

 6883 CONTINUE

C-----DE NUMERACION DARIO A LA MIA (COORDS)
      CALL DADO(R,NT,NATSOL,IZ)
C-----Pone las cargas en el solvente(en doda se cambiaron)
C-----DEBE SER IZZ=1

C=======ELFIELD=1 ===========================
      IF(ELFIELD.EQ.1) THEN
      DO I=1+NATOM,NWAT+NATOM
      PC(I)  = PC0(I)
      ENDDO
      GOTO 6884
      ENDIF
ccccccccccccccccccccccccccccccccccccccccccccccccc

C      DO I=1+NATOM,NWAT+NATOM
C      I1= I+NWAT
C      I2= I+2*NWAT
 
C      PC(I)  = PC0(I)
C      PC(I1) = PC0(I1)
C      PC(I2) = PC0(I2)

C      ENDDO

 6884 CONTINUE

      ELSEIF (NPAS.GT.0.AND.MOD(ITEL,NPAS).NE.0) THEN
      
    
C-- Si no hace en CADA paso un calculo de estructura electronica:
C-- pero SI calcula en CADA paso integral de RHO*V, Zj*Pc(i)/|Rj-Ri| 
C-- y potencial LJ entre Ri y Rj:

      CALL DODA(NATSOL,R,NT,IZ)
      
      CALL INTSOL(NORM,natom,Nsol,natsol,r,Nuc,Iz,M,Md,ncont,nshell,
     > c,a,pc,RMM,E1s,FQQ,IT,ITEL,NIN,IPR1,EAC,NPAS)
      

      DO I=1,NPART
       DO J=1,3
       f1(I,J)=0.0D0
      ENDDO
      ENDDO

      CALL INTSOLG(NORM,natom,Nsol,natsol,r,Nuc,Iz,M,Md,
     > ncont,nshell,c,a,pc,RMM,f1,FXH11,FYH11,FZH11)
      FXH11 =  f1(1,1) - FXH2
      FYH11 =  f1(1,2) - FYH2
      FZH11 =  f1(1,3) - FZH2
       
C-----CAMBIO UNIDADES FZA. SOBRE SOLUTO
      DO I=1,NATOM
      FX(I)=FX(I)-F1(I,1)/A0*HH
      FY(I)=FY(I)-F1(I,2)/A0*HH
      FZ(I)=FZ(I)-F1(I,3)/A0*HH
c      write(*,*)'f1 ',i,(f1(i,j)/a0*hh,j=1,3)
      ENDDO

C-----ACOPLAMIENTO COULOMBIANO SOLUTO-SOLVENTE POR ATOMO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(ELFIELD.EQ.1) THEN
      I1=NATOM
      DO I=1+NATOM,NPART
      I1=I1+1   !   O
      VCSS(I1)=VCSS(I1)*FACTOR + EAC(I)*HH*FACTOR
      VCSSO = VCSSO + VCSS(I1)
      ENDDO
      GOTO 6885
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      I1=NATOM
      DO I=1+NATOM,NPART,3
      I1=I1+1   !   O
      VCSS(I1)=VCSS(I1)*FACTOR + EAC(I)*HH*FACTOR
      VCSSO = VCSSO + VCSS(I1)
      ENDDO

      DO I=1+NATOM,NPART,3
      I1=I1+1   !   H
      VCSS(I1)=VCSS(I1)*FACTOR + EAC(I+1)*HH*FACTOR
      VCSSH = VCSSH + VCSS(I1)
      ENDDO

      DO I=1+NATOM,NPART,3
      I1=I1+1   !   H
      VCSS(I1)=VCSS(I1)*FACTOR + EAC(I+2)*HH*FACTOR
      VCSSH = VCSSH + VCSS(I1)
      ENDDO

 6885 CONTINUE

      VCSST = VCSST + VCSSO + VCSSH

C-----TRANSFORMA NUMERACION Y UNIDADES FZA. SOBRE SOLVENTE 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(ELFIELD.EQ.1) THEN
      NN=NATOM
      DO N=NATOM+1,NPART            
      NN=NN+1
      FX(NN)=FX(NN)-F1(N,1)/A0*HH
      FY(NN)=FY(NN)-F1(N,2)/A0*HH
      FZ(NN)=FZ(NN)-F1(N,3)/A0*HH
      FQ(NN)=FQ(NN)+FQQ(N)/EE*HH
      ENDDO
      GOTO 6886
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      NN=NATOM
      DO N=NATOM+1,NATOM+NWAT*NATSOL,3
      NN=NN+1
      FX(NN)=FX(NN)-F1(N,1)/A0*HH
      FY(NN)=FY(NN)-F1(N,2)/A0*HH
      FZ(NN)=FZ(NN)-F1(N,3)/A0*HH
      FQ(NN)=FQ(NN)+FQQ(N)/EE*HH
      ENDDO
 
      DO N=NATOM+1,NATOM+NWAT*NATSOL,3
      NN=NN+1
      FX(NN)=FX(NN)-F1(N+1,1)/A0*HH
      FY(NN)=FY(NN)-F1(N+1,2)/A0*HH
      FZ(NN)=FZ(NN)-F1(N+1,3)/A0*HH
      FQ(NN)=FQ(NN)+FQQ(N+1)/EE*HH
      ENDDO
 
      DO N=NATOM+1,NATOM+NWAT*NATSOL,3
      NN=NN+1
      FX(NN)=FX(NN)-F1(N+2,1)/A0*HH
      FY(NN)=FY(NN)-F1(N+2,2)/A0*HH
      FZ(NN)=FZ(NN)-F1(N+2,3)/A0*HH
      FQ(NN)=FQ(NN)+FQQ(N+2)/EE*HH
      ENDDO

 6886 CONTINUE

      CALL DADO(R,NT,NATSOL,IZ)

      ENDIF
9000  CONTINUE
      IF (IEWLD.EQ.1) CALL FURIM
      Vcerca = Vcerca*FACTOR
      Vlejos = Vlejos*FACTOR
      VCT  = VCT * FACTOR
      VLJ  = VLJ * FACTOR
      VOO  = VOO * FACTOR 
      VCOO = VCOO* FACTOR 
      VOH  = VOH * FACTOR 
      VHH  = VHH * FACTOR 
      VTOT = VTOT* FACTOR
      VSELF= VSELF*FACTOR
      VCQC = VCQC* FACTOR
      VLJQC= VLJQC*FACTOR

C------------------------------------------------------------------C
      FXT = ZERO
      FYT = ZERO
      FZT = ZERO

 
      DO 201 I = 1, NPART
      FXT = FXT + FX (I)
      FYT = FYT + FY (I)
      FZT = FZT + FZ (I)
201   CONTINUE
c      write(*,*)'Fza total: ',abs(fxt*fxt+fyt*fyt+fzt*fzt)
c      write(*,*)'           ',fxt,fyt,fzt

      IF (DSQRT(FXT*FXT+FYT*FYT+FZT*FZT).GT.1.D-04)THEN
      WRITE (6,*) ' FZA TOTAL NE ZERO EN MAIN  ', FXT, FYT,FZT
      ENDIF
C------------------------------------------------------------------C


C-----MUEVE COORDENADAS:
C-----DE TODOS LOS ATOMOS SI LLAMO A DFT,
C-----SOLO EL SOLVENTE SI NO LLAMO A DFT (SOLUTO CONGELADO)

        IF(BER1.NE.BER2) THEN
        	IF(BER1.LT.BER2)THEN
        	BERDK = (BER2-BER1)/BERSTEP
        	KBER = KBER + BERDK
        	ELSE
        	BERDK = (BER1-BER2)/BERSTEP
        	KBER = KBER - BERDK
        	ENDIF
        ELSE
        KBER=BER2
        ENDIF

      QSEK = ZERO
      QSEP = ZERO

      IF(ICMOT.EQ.1)THEN
         IF(NDFT*NPAS.GT.0.AND.MOD(ITEL,NPAS).EQ.0)THEN
           CALL CVTMOT(NATSOL)
         ELSEIF(NDFT.NE.1.OR.MOD(ITEL,NPAS).NE.0)THEN
           CALL CVTMOT1(NATSOL)
         ENDIF
                    QSEP = GDF*TEMPRQ*S0*FACTOR
                    QSEK = QS*SD0*SD0
 
      ELSEIF(ICMOT.EQ.0)THEN
            IF(NDFT*NPAS.GT.0.AND.MOD(ITEL,NPAS).EQ.0)THEN
              CALL CVEMOT(NATSOL)
            ELSEIF(NDFT.NE.1.OR.MOD(ITEL,NPAS).NE.0)THEN
              CALL CVEMOT1(NATSOL)
            ENDIF

      ELSEIF(ICMOT.EQ.2)THEN
C-----SI ICMOT.EQ.2 HACE STEEPEST DESCENT CON TODOS LOS NUCLEOS
            CALL OPTIM(NATSOL)

      ELSEIF(ICMOT.EQ.3)THEN
            IF(NDFT*NPAS.GT.0.AND.MOD(ITEL,NPAS).EQ.0)THEN
              CALL CVEMOTN(NATSOL)
            ELSEIF(NDFT.NE.1.OR.MOD(ITEL,NPAS).NE.0)THEN
              CALL CVEMOTN1(NATSOL)
            ENDIF
 
       ELSEIF(ICMOT.EQ.4)THEN
C            IF(NDFT*NPAS.GT.0.AND.MOD(ITEL,NPAS).EQ.0)THEN
               CALL CVEMOTNN(NATSOL)
C            ELSEIF(NDFT.NE.1.OR.MOD(ITEL,NPAS).NE.0)THEN
C              CALL CVEMOTN1(NATSOL)
 
       ELSEIF(ICMOT.EQ.5)THEN
c              write(*,*) 'ANTES',X(6),Y(6),Z(6)
C               CALL CVEMOTNNN(NATSOL)
c              write(*,*) 'DESPUES',X(6),Y(6),Z(6)
           CONTINUE
       ENDIF

C-----MUEVE CARGAS SIST.CLASICO
      IF(NWAT.EQ.0)GOTO 7447
      IF(IFLUC.NE.1)GOTO 3234

      QSEKQ = ZERO
      QSEPQ = ZERO

      IF (IQMOT.EQ.1)THEN
CQTMOTB ==> a temperatura cte con Berendsen/QSQ=>masa termostato

                    CALL QTMOTB

      ELSEIF(IQMOT.EQ.0)THEN

                  CALL QEMOT

      ELSEIF(IQMOT.EQ.2)THEN
C-----SI IQMOT.EQ.2 HACE STEEPEST DESCENT CON LAS CARGAS PARCIALES
                  CALL OPTIMQ

      ELSE

      CONTINUE
                 
      ENDIF

3234  CONTINUE
7447  CONTINUE


      TEMPSQ= QSEKQ/BOLTZF
      TEMPS = QSEK/BOLTZF
      QSEKQ = PTFIVE*QSEKQ*FACTOR
      QSEK = PTFIVE*QSEK*FACTOR


C-------------------------------------------------------C
C-----SI ESTA OPTIMIZANDO: CONTROLA LAS FUERZAS
      FQT2 = ZERO
      FXT2 = ZERO
      FYT2 = ZERO
      FZT2 = ZERO
      FQT2 = ZERO

      DO 211 I = 1, NPART
      FXT2 = FXT2 + FX(I)*FX(I)
      FYT2 = FYT2 + FY(I)*FY(I)
      FZT2 = FZT2 + FZ(I)*FZ(I)
      FQT2 = FQT2 + FQ(I)*FQ(I)
211   CONTINUE
      
      EMOD = (FXT2 + FYT2 + FZT2)
      EMODQ = (FQT2)
      EMOD = DSQRT(EMOD)*FACTOR*A0
      EMODQ = DSQRT(EMODQ)*FACTOR*EE
C--------------------------------------------------------C


C-----TERMALIZA NUCLEOS
      IF (NEWV.GT.0.AND. MOD(ITEL,NEWV).EQ.0) THEN

          IF(NSCAL.EQ.0.AND.ITEL.EQ.0) GOTO 353
                
		IF(NDFT*NPAS.EQ.1) CALL TEPCHR(NATSOL)
		IF(NDFT*NPAS.NE.1) CALL TEPCHR2(NATSOL)
                TEMPS = ZERO

      ENDIF
353   CONTINUE


C-----TERMALIZA CARGAS SIST.CLASICO
      IF (NEWQ.GT.0.AND. MOD(ITEL,NEWQ).EQ.0) THEN
               
                CALL TEPQ
                TEMPSQ = ZERO

      ENDIF

C-----EJES PPALES DE INERCIA Y MOM. DIPOLAR TOTAL
C          IF(NDIP.GT.0)THEN
C     CALL UNGLES(NATSOL,IT,NIN,NDIP)      
C     CALL DIPOLOS(IT,NIN,NDIP,DV1,DV2,DV3,
C    >  ux,uy,uz)
C          ENDIF

551   CONTINUE

      RT = ONE / DBLE (IT)
      VOO = VLJ + VCOO
      EKIN = PTFIVE * TEMPAV * FACTOR
c      WRITE(*,*) 'EKIN,TEMPAV,FACTOR',EKIN,TEMPAV,FACTOR
      EKINQ =PTFIVE * TEMPOL * FACTOR
      IF(NWAT.EQ.0)THEN
      EK = EKIN
      ELSE
      EK = EKIN + EKINQ
      ENDIF
      ECSELF = VSELF - DBLE(NWAT)*EGP*FACTOR
      EINTQC = VCQC + VLJQC
      IF(IFLUC.EQ.1)THEN
      EPOT = VCT + VLJ + ENEFUR*FACTOR + ECSELF
      ELSE  
      EPOT = VCT + VLJ + ENEFUR*FACTOR 
      ENDIF
      EPOT1= EPOT + EINTQC 
      POT = (EKS+E1s)*HH*FACTOR + EPOT1
      ETOTL = EK + POT
      TE = ETOTL+(QSEK+QSEP+QSEKQ+QSEPQ)
c      WRITE(*,*)'TEMPAV,GDF',TEMPAV,GDF
      TEMPAV = TEMPAV/GDF
c     WRITE(*,*) 'TEMPAV',TEMPAV
       if(GDFSLT.eq.zero)then
      TEMPSLT=TEMPSLT
      TEMPSLV=TEMPSLV/gdf
       else
      TEMPSLV=TEMPSLV/GDFSLV
      TEMPSLT=TEMPSLT/GDFSLT
       endif
      TEMPOL = TEMPOL/GDQ

      
      AC (1)  = AC(1)  + TEMPOL
      AC (2)  = AC(2)  + TEMPOL*TEMPOL 
      AC (3)  = AC(3)  + TEMPAV
      AC (4)  = AC(4)  + TEMPAV* TEMPAV
      AC (7)  = AC(7)  + S0
      AC (8)  = AC(8)  + QSEP
      AC (9)  = AC(9)  + TEMPS
      AC (10) = AC(10) + ETOTL
      AC (11) = AC(11) + ETOTL * ETOTL
      AC (12) = AC(12) + TE
      AC (13) = AC(13) + TE * TE
      AC (14) = AC(14) + EPOT
      AC (15) = AC(15) + EPOT **2
      AC (16) = AC(16) + VLJ
      AC (17) = AC(17) + VLJ * VLJ
      AC (18) = AC(18) + VCOO
      AC (19) = AC(19) + VCOO * VCOO
      AC (20) = AC(20) + VOO
      AC (21) = AC(21) + VOO * VOO
      AC (22) = AC(22) + VOH
      AC (23) = AC(23) + VOH * VOH
      AC (24) = AC(24) + VHH
      AC (25) = AC(25) + VHH * VHH
      AC (26) = AC(26) + VCT
      AC (27) = AC(27) + VCT * VCT
      AC (28) = AC(28) + VSELF
      AC (29) = AC(29) + VSELF * VSELF
      AC (30) = AC(30) + ECSELF
      AC (31) = AC(31) + ECSELF * ECSELF
      AC (32) = AC(32) + SUMFX
      AC (33) = AC(33) + SUMFY
      AC (34) = AC(34) + SUMFZ
      AC (35) = AC(35) + SUMF
      AC (36) = AC(36) + EINTQC
      AC (37) = AC(37) + EINTQC * EINTQC
      AC (38) = AC(38) + VCQC
      AC (39) = AC(39) + VCQC*VCQC
      AC (40) = AC(40) + VLJQC
      AC (41) = AC(41) + VLJQC*VLJQC
      AC (44) = AC(44) + EPOT1
      AC (45) = AC(45) + EPOT1 **2
      AC (46) = AC(46) + POT
      AC (47) = AC(47) + POT*POT 
      AC (48) = AC(48) + DV1
      AC (49) = AC(49) + DV1*DV1
      AC (50) = AC(50) + DV2
      AC (51) = AC(51) + DV2*DV2
      AC (52) = AC(52) + DV3
      AC (53) = AC(53) + DV3*DV3
      AC (54) = AC(54) + DIPT
      AC (55) = AC(55) + DIPT*DIPT
      AC (56) = AC(56) + VCSSO
      AC (57) = AC(57) + VCSSO*VCSSO
      AC (58) = AC(58) + VCSSH
      AC (59) = AC(59) + VCSSH*VCSSH
      AC (60) = AC(60) + VCSST
      AC (61) = AC(61) + VCSST*VCSST
      AC (62) = AC(62) + FZ5
      AC (63) = AC(63) + FZ5*FZ5
      AC (64) = AC(64) + FZ1
      AC (65) = AC(65) + FZ1*FZ1
      AC (66) = AC(66) + FZ2
      AC (67) = AC(67) + FZ2*FZ2 
      AC (68) = AC(68) + Vcerca 
      AC (69) = AC(69) + Vlejos 


C---- ESTADISTICA PARA LAS POBLACIONES DE MULLIKEN (Q(I))
C---- Y PARA LAS CARGAS FLUCTUANTES (PC(I))   
c----------------c
      NN = 67
c----------------c
c      DO 22 I=1,NPART
c      IF(I.LE.NATOM)THEN
c      AC(NN+I) = AC(NN+I) + Q(I)
c      AC(I+NN+NPART) = AC(I+NN+NPART) + Q(I)*Q(I)
c      ELSE
c      AC(NN+I) = AC(NN+I) + PC(I)/EE
c      AC(NN+I+NPART) = AC(NN+I+NPART) + PC(I)/EE*PC(I)/EE
c      ENDIF
c22    CONTINUE   


      NREF=NATOM+NWAT+NN+2*NPART
c      IF(NREF.GT.500)THEN
c      WRITE(*,*)'DANGER! DIMENSION AC(I).GT.500 ',NREF
c      WRITE(*,*)'CUIDADO CON LOS ACUMULADORES DE CARGAS Y ETC'
c      STOP
c      ENDIF

      DO 290 I = 1, 500
      AV(I) = AC(I) * RT
290   CONTINUE
      
      IF (MOD((IT-NIN),IPR).EQ.0)THEN
      
          IF(MOD(ITEL,NPAS).NE.0.OR.NDFT.NE.1)THEN

             IF(ICMOT.NE.2.AND.IQMOT.NE.2)THEN
      WRITE (6,350)ITEL,TE,EK,POT,E1s*HH*FACTOR,EKS*HH*FACTOR,
     & TEMPAV,TEMPOL
             ENDIF
          ELSE

      WRITE(6,*)
      WRITE(6,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      WRITE(6,*)'STEP No:  ',ITEL
      WRITE(6,*)'TEMPERATURA NUCLEOS           ',TEMPAV
      WRITE(6,*)'TEMPERATURA CARGAS            ',TEMPOL
      WRITE(6,*)'ENERGIA POLARIZACION          ',ECSELF
      WRITE(6,*)'ENERGIA LJ SOLVENTE           ',VLJ
      WRITE(6,*)'ENERGIA ES SOLVENTE           ',VCT
      WRITE(6,*)'ENERGIA POT. SOLVENTE CLASICO ',EPOT
      WRITE(6,*)'ENERGIA VLJQC                 ',VLJQC
      WRITE(6,*)'ENERGIA VCQC                  ',VCQC
      WRITE(6,*)'ENERGIA POT. VLJQC + VCQC     ',EINTQC
      WRITE(6,*)'ENERGIA POT. CLAS+ INTQC      ',EPOT1
      WRITE(6,*)'ENERGIA RHO*V                 ',E1s
c      WRITE(6,*)'ENERGIA RHO*V                 ',E1s*HH*FACTOR
c      WRITE(6,*)'ENERGIA KS                    ',EKS*HH*FACTOR
      WRITE(6,*)'ENERGIA KS                    ',EKS
      WRITE(6,*)'ENERGIA POTENCIAL TOTAL       ',POT
      WRITE(6,*)'ENERGIA CINETICA              ',EK
      WRITE(6,*)'HAMILTONIANO                  ',TE
      WRITE(6,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      WRITE(6,*)

         ENDIF
     
      ENDIF

c                     *******
                      NPR=10 
c                     *******
           IF(MOD((IT-NIN),NPR).EQ.0)THEN

C-----ENERGIAS Y TEMPERATURAS EN F. DEL TIEMPO
      IF(ICMOT.NE.2)THEN
      WRITE(19,86)ITEL,TE,EK,POT,EINTQC,ECSELF,TEMPAV,TEMPOL,E1s*HH*FACTOR,
     & EKS*HH*FACTOR
      WRITE(80,'(i6,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3)')ITEL,TEMPAV,
     & TEMPOL,TEMPSLT,TEMPSLV
      CALL FLUSH(80)
c      WRITE(81,*)ITEL,DIST,DIST1,DIST2
      WRITE(82,*)ITEL,E1s*HH*FACTOR,VCQC
      WRITE(83,*)ITEL,(VCSS(I),I=1+NATOM,NATOM+NWAT)
      WRITE(84,*)ITEL,(VCSS(I),I=1+NATOM+NWAT,NPART)
      WRITE(85,*)ITEL,VCSSO,VCSSH,VCSST
      ELSEIF(ICMOT.EQ.2)THEN
C-----SI OPTIMIZA Q Y R ESCRIBE ENERGIA TOTAL, POTENCIAL Y GRADIENTES
      WRITE(19,86)ITEL,TE,EMOD,EMODQ
      IF(NSPECQ.EQ.0)THEN
      WRITE(6,86)ITEL,TE,EMOD,EMODQ
      ENDIF
      ENDIF

           ENDIF
     
C-----IMPRIME ARCHIVO FORT.15 CADA "IP15" PASOS:
c     write(*,*)'ANTES DE 15',X(6),Y(6),Z(6)
 
         IF(MOD((IT-NIN),IP15).EQ.0)THEN
      WRITE(15,*)NPART,ITEL
      WRITE(15,*)
      DO I=1,NPART
        IF(NDFT.NE.1)THEN
      WRITE(15,101)AT(I),X(I),Y(I),Z(I),PC(I)/EE
        ELSE
      IF(I.LE.NATOM)THEN
      WRITE(15,101)AT(I),X(I),Y(I),Z(I),Q(I)
      ELSE
      WRITE(15,101)AT(I),X(I),Y(I),Z(I),PC(I)/EE
      ENDIF
        ENDIF
      ENDDO
         ENDIF

c     write(*,*)'DESPUES DE 15',X(6),Y(6),Z(6)

C --  ESCRIBE EN CADA PASO: COORDS. Y VELOCIDADES SOLUTO
      WRITE(18,*)NATOM,ITEL
      WRITE(20,*)ITEL
      WRITE(18,*)
      DO I=1,NATOM
      IF(NDFT.NE.1)WRITE(18,101)AT(I),X(I),Y(I),
     & Z(I),PC(I)/EE
      IF(NDFT.EQ.1)WRITE(18,101)AT(I),X(I),Y(I),
     & Z(I),Q(I)
      WRITE(20,101)AT(I),VX(I),VY(I),VZ(I)
      ENDDO
C-----END FORT.15 y FORT.18 y FORT.20

      IF (MOD(ITEL,NSAVE).EQ.0.AND.NSAVE.GT.0)THEN
      WRITE (8,*)IDUM,X99
      WRITE (8,*)ITEL
      WRITE (8,*)X,Y,Z,X0,Y0,Z0,VX,VY,VZ,VX0,VY0,VZ0
      WRITE (8,*)ITEL
      WRITE (8,*)ST,S0,SD,SD0
      WRITE (8,*)STQ,S0Q,SDQ,SD0Q
      WRITE (8,*)ITEL
      WRITE (8,*)PC,PC0,VQ
      WRITE (8,*)ITEL
      WRITE (8,*)KG,KG1
      WRITE (8,*)HISTO
      WRITE (8,*)HISTO2
      WRITE (8,*)QK
      WRITE (8,*)AC
      REWIND 8

      ENDIF

C-----LLAMA A 'GRILLA' PARA G(R) ELECTRONICO(SOLV)
c      CALL GRILLA(ntq,Q)

390   CONTINUE


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
CC         FINAL DEL LOOP
CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C-----CALCULA E IMPRIME EN FORT.16 Y FORT.17 G(R)
      IF(IRDF.EQ.1)THEN
      CALL PRINTER(RT,HISTO,HISTO2,DELP,DELPZ)
      ENDIF

C-----IMPRIME FINAL EN FILE8.OUT        
      WRITE (8,*)IDUM,X99
      WRITE (8,*)ITEL
      WRITE (8,*)X,Y,Z,X0,Y0,Z0,VX,VY,VZ,VX0,VY0,VZ0
      WRITE (8,*)ITEL
      WRITE (8,*)ST,S0,SD,SD0
      WRITE (8,*)STQ,S0Q,SDQ,SD0Q
      WRITE (8,*)ITEL
      WRITE (8,*)PC,PC0,VQ
      WRITE (8,*)ITEL
      WRITE (8,*)KG,KG1,HISTO
      WRITE (8,*)HISTO2
      WRITE (8,*)QK
      WRITE (8,*)AC

C----IMPRIME LOS INPUT=FILE8.out PARA JARZINSKY
c         
c       IF(ITEL.EQ.IPINPUT0) THEN
c
cc        open(unit=IPINPUT0,file=file8.IPINPUT0)
c
c         IPINPUT0 = IPINPUT0+IPINPUT
c
c         WRITE (,*)IDUM,X99
c         WRITE (,*)ITEL
c         WRITE (,*)X,Y,Z,X0,Y0,Z0,VX,VY,VZ,VX0,VY0,VZ0
c         WRITE (,*)ITEL
c         WRITE (,*)ST,S0,SD,SD0
c         WRITE (,*)STQ,S0Q,SDQ,SD0Q
c         WRITE (,*)ITEL
c         WRITE (,*)PC,PC0,VQ
c         WRITE (,*)ITEL
c         WRITE (,*)KG,KG1,HISTO
c         WRITE (,*)HISTO2
c         WRITE (,*)QK
c         WRITE (,*)AC
c
c       ENDIF

C-----REESCRIBE LA FOTO (tt.alc) "FINAL" 
        WRITE(41,50) NPART
        WRITE(41,*) 
        DO I=1,NPART
        IF(NDFT.EQ.1)THEN
         IF(I.LE.NATOM)THEN
          WRITE(41,101) AT(I),X(I),Y(I),Z(I),Q(I)
         ELSE
          WRITE(41,101) AT(I),X(I),Y(I),Z(I),PC(I)/EE
         ENDIF
        ELSE
        WRITE(41,101) AT(I),X(I),Y(I),Z(I),PC(I)/EE
        ENDIF
        ENDDO
        WRITE(41,*)'ITEL: ',ITEL
C-----FIN FOTO  

      WRITE (6,*)

      IF(IUNID.EQ.1)WRITE(6,*)'UNIDADES: ADIMENSIONALES, FACTOR DE 
     & NORMALIZACION= 1/NKT'
      IF(IUNID.EQ.2)WRITE(6,*)'UNIDADES: KELVIN'
      IF(IUNID.EQ.3)WRITE(6,*)'UNIDADES: Kcal/mol'
      IF(IUNID.EQ.4)WRITE(6,*)'UNIDADES: HARTREES'

      IF(ICMOT.EQ.2)GOTO 909
     
      WRITE (6,*)
      WRITE (6,*) ' .... AFTER  ',ITEL,' STEPS HAVE BEEN ACCUMULATED'
      WRITE (6,*)
      WRITE (6,*)'<HAMILTONIAN>', AV(12),' <DH>',DSQRT(AV(13)-AV(12)**2)
      WRITE (6,*)
      WRITE (6,*) '<TEMPERATURE>',AV(3),' <DT>',DSQRT(AV(4)-AV(3)**2)
      WRITE (6,*)
      WRITE (6,*) '<TEMPOL> ',AV(1),' <DT>',DSQRT(AV(2)-AV(1)**2)  
      WRITE(6,*)
      WRITE (6,*) '<POT CLAS >',AV(14),' <DP>',DSQRT(AV(15)-AV(14)**2)
      WRITE (6,*)
      WRITE (6,*) '<POT L-J >',AV(16),' <DP>',DSQRT(AV(17)-AV(16)**2)
      WRITE (6,*)
      WRITE (6,*) '<POT COUL >',AV(26),' <DV>',DSQRT(AV(27)-AV(26)**2)
      WRITE (6,*)
c      WRITE (6,*) '<POT VCQC>',AV(38),'  <VCQC**2>',AV(39),
c     &  '  <D VCQC 2>',DSQRT(AV(39)-AV(38)**2)
c      WRITE (6,*)
c      WRITE (6,*) '<POT LJQC>',AV(40),'  <VLJQC**2>',AV(41),
c     &  '  <D VLJQC 2>',DSQRT(AV(41)-AV(40)**2)
c      WRITE (6,*)
c      WRITE (6,*) '<POT INTQC >',AV(36),'  <EIQC**2>',AV(37),
c     &  '  <D  2>',DSQRT(AV(37)-AV(36)**2)
c      WRITE (6,*)
      WRITE (6,*) '<POT (C)+(QC) >',AV(44),' <DP>',
     & DSQRT(AV(45)-AV(44)**2)
      WRITE (6,*)
      WRITE (6,*) '<POT TOTAL >',AV(46),' <DP>',DSQRT(AV(47)-AV(46)**2)
      WRITE (6,*)
c      WRITE (6,*) '< DV1 >',AV(48),'  <D DV1 >',DSQRT(AV(49)-AV(48)**2)
c      WRITE (6,*)
c      WRITE (6,*) '< DV2 >',AV(50),'  <D DV2 >',DSQRT(AV(51)-AV(50)**2)
c      WRITE (6,*)
c      WRITE (6,*) '< DV3 >',AV(52),'  <D DV3 >',DSQRT(AV(53)-AV(52)**2)
c      WRITE (6,*)
      WRITE (6,*) '< V >',AV(54),'  <D V >',DSQRT(AV(55)-AV(54)**2)
      WRITE (6,*) 
      WRITE (6,*) '<VCSSO>',AV(56),'   <Delta>',DSQRT(AV(57)-AV(56)**2)
      WRITE (6,*)
      WRITE (6,*) '<VCSSH>',AV(58),'   <Delta>',DSQRT(AV(59)-AV(58)**2)
      WRITE (6,*)
      WRITE (6,*) '<VCSST>',AV(60),'   <Delta>',DSQRT(AV(61)-AV(60)**2)
      WRITE (6,*)

c      WRITE (6,*) '<FZ>',AV(62),'   <Delta>',DSQRT(AV(63)-AV(62)**2)
c      WRITE (6,*)
c      WRITE (6,*) '<F1>',AV(64),'   <Delta>',DSQRT(AV(65)-AV(64)**2)
c      WRITE (6,*)
c      WRITE (6,*) '<F2>',AV(66),'   <Delta>',DSQRT(AV(67)-AV(66)**2)
c      WRITE (6,*)
       write (6,*) 'nano capo, Vcerca y Vlejos', AV(68), AV(69)


c      OPEN(7,FILE='mull.out')
c      WRITE (7,*)
c      WRITE(7,*)'<MULLIKEN POP Y CARGAS: >'
c      DO I=1,NPART
c      WRITE(7,612)I,IZ(I),AV(NN+I),DSQRT(AV(NN+I+NPART)-AV(NN+I)**2)
c      ENDDO

909   CONTINUE


      WRITE(*,*)
      WRITE(*,*)'----- FIN -----'
      date='date'
      CALL SYSTEM(date)



86    FORMAT(I10,4X,10G20.9)
87    FORMAT(2X,10G15.7)
45    FORMAT(2X,I7,10F12.6)
778   FORMAT(2X,I8,10G18.7)
779   FORMAT(2X,I8,20G18.7)
223   FORMAT(2X,I5,8X,4G15.7)   
612   FORMAT(3X,2I5,5X,4G18.7)
250   FORMAT(4X,'ITER',4X,'ENERGY1',6X,'EKINQ',6X,
     &'EPOT1',3X,' VLJ  ',5X,'VCOUL ',4X,'ESELF',6X,'EINTQC',4X,
     &'TEMPAV',4X,'TEMPOL',/
     &,'--------------------------------------------------------')
466   FORMAT(I7,10F9.2)
467   FORMAT(I7,20F9.2)
350   FORMAT (I7,2X,11G12.5)
91    FORMAT (8I10)
92    FORMAT (4G19.7)
93    FORMAT (4G19.10)
78    FORMAT(2X,'D (H Br) ',2I4,2G18.9)
77    FORMAT(2X,'D H)(O   ',2I4,2G18.9)
79    FORMAT(2X,'D (H O)  ',2I4,2G18.9)
50    FORMAT(2X,I3,1X,'ATOMS')
100   FORMAT(2X,I3,1X,A5,4F9.4)
101   FORMAT(5X,A5,4F9.4)
102   FORMAT(2X,A5,4F15.6)



      END


