      SUBROUTINE INICIO(NATSOL,NDIP,IDIPCOR,PMAX,PZMAX)
      INCLUDE 'COMM'
      
      INTEGER SPC  
      COMMON /tipsol/SPC
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C       ICON: DEFINE EL TIPO DE CORRIDA, 0: COMIENZA CON ORDEN
C             1: LEE DE 8 Y CONTINUA EL AVERAGE
C            -1: LEE SOLO LA CONFIGURACION DE 8 Y ANULA LOS AVERAGES
C       NTIME: NUMERO DE PASOS DE LA SIMULACION
C       NSCAL: DECIDE SI A t = 0 (INICIO) RANDOMIZA VELOCIDADES
C       NSAVE: CADA TANTO REWIND Y REESCRIBE EN EL REGISTRO 8
C       NDIP: CALCULA O NO MOMENTOS DIPOLARES (=1 =>SI )
C       NANG: CONTROLA INPRESION EN ARCHIVO 11
C       IRDF: CALCULA O NO g(r), 1:CALCULA G(r), SINO NO
C       ICMOT: DEFINE EL TIPO DE DINAMICA: MICROCANONICA O CANONICA
C       IQMOT: DEFINE DINAMICA DE LAS CARGAS, MICRO O CANONICA 
C       NEWV: CADA CUANTOS PASOS VUELVE A RANDOMIZAR LAS TEMPERATURAS
C       NEWQ: CADA CUANTOS PASOS RANDOMIZA TEMP. DE LAS CARGAS
C       TEMPRQ: TEMPERATURA DEL TERMOSTATO CONFIGURACIONAL
C       TEMPRQQ: TEMP.TERMOSTATO CARGAS FLUCTUANTES
C       IPR: CADA CUANTOS PASOS IMPRIME
C       IEWLD: EWALD SUM PARA LAS FUERZAS (FSPHERCUT-FURIM)
C              SI  IEWLD=0, NO LLAMA FURIM 
C       XFAX: DEFINE EL RCUTOFF(RCT) 
C       ICLSTR: DEFINE SI ES UN CLUSTER O FASE CONDENSADA
C                1,CLUSTER     SINO ,BULK 
C       UMBRE: define si hace Umbrella sampling
C       CKUM y RUM: Constante y distancia del umbrella
C       US1 y US2: Atomos cuanticos que usare en el umbrella 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      OPEN (51,FILE='constraint.dat')
      OPEN (8,FILE='file8.out')
      OPEN (12,FILE='ini.xyz')


      READ (51,*)
      READ (51,*) ICON, NTIME, NSCAL,NSAVE,NDIP,IRDF,ICMOT,NEWV,
     &  IPR,IP15,UMBRE,BER1,BERSTEP,BER2  
      READ (51,*)
      READ (51,*) IDUM,NEWQ,IQMOT,IZZ,IFLUC,IPR1,IS0,IS0Q,
     & CKUM,RUM,JVEL,US1,US2,US3,US4
       write(*,*) 'umbrella', CKUM, RUM, jvel,US1,US2,US3,US4

c     READ (51,*) 
c     READ (51,*) ITERM,IPINPUT
      WRITE (6,55)ICON,NTIME,NSCAL,NSAVE
      WRITE (6,59)NDIP,IRDF,ICMOT,NEWV,IPR,NANG
      READ (51,*)
      READ (51,*) NSPECQ,NDFT,NPAS
      READ (51,*) 
      READ (51,*) (NNAT(I),I=1,NSPECQ)

c--ojo, soluto totalmente congelado:
      IF(NDFT.NE.1)THEN
      READ (51,*)
      READ (51,*)(PC1(I),I=1,NSPECQ)
      write(*,*)' >>>>>SOLUTO TOTALMENTE CONGELADO!!<<<<<<'
      ELSE
c---sino:
      READ(51,*)
      READ(51,*)
      ENDIF
c---
      READ (51,*) 
      READ (51,*) (WWM(I),I=1,NSPECQ)
      TOLE=0.2D0
      DO I=1,NSPECQ
      IF(ABS(WWM(I)-16.D0).LT.TOLE)AT(I)='O    '
      IF(ABS(WWM(I)-35.D0).LT.TOLE)AT(I)='Cl   '
      IF(ABS(WWM(I)-1.D00).LT.TOLE)AT(I)='H    '
      IF(ABS(WWM(I)-2.D00).LT.TOLE)AT(I)='D    '
      IF(ABS(WWM(I)-14.D0).LT.TOLE)AT(I)='N    '
      IF(ABS(WWM(I)-39.1D0).LT.TOLE)AT(I)='K    '
      IF(ABS(WWM(I)-12.D0).LT.TOLE)AT(I)='C    '
      IF(ABS(WWM(I)-32.D0).LT.TOLE)AT(I)='S    '
      IF(ABS(WWM(I)-126.9D0).LT.TOLE)AT(I)='I    '
      IF(ABS(WWM(I)-6.941D0).LT.TOLE)AT(I)='Li   '
      IF(ABS(WWM(I)-26.98D0).LT.TOLE)AT(I)='Al   '
      ENDDO
      WRITE(*,1000)(AT(I),I=1,NSPECQ)
      READ (51,*)
      READ (51,*) NATSOL
      WRITE (6,70)NATSOL
      READ (51,*)
      READ (51,*) (NNAT(I),I=NSPECQ+1,NSPECQ+NATSOL)
      WRITE (6,656) (NNAT(I),I=NSPECQ+1,NSPECQ+NATSOL)
      READ (51,*)
      READ (51,*) (WWM(I),I=NSPECQ+1,NSPECQ+NATSOL+1)
c      WRITE (6,56) (WWM(I),I=1,NSPECQ+NATSOL+1)
      READ (51,*)
      READ (51,*) TEMPRQ,TEMPRSV,TEMPRQQ,DELTAT,QS, QSQ
      WRITE (6,57) TEMPRQ,TEMPRQQ,DELTAT
      WRITE(6,60) QS,QSQ
      READ (51,*)
      READ(51,*)  (SIGMA(I),I=NSPECQ+1,NSPECQ+NATSOL)
      WRITE (6,111)SIGMA(NSPECQ+1)
      READ (51,*)
      READ(51,*)  (EPS(I),I=NSPECQ+1,NSPECQ+NATSOL)
      WRITE (6,112)EPS(NSPECQ+1) 
      READ (51,*)
      READ(51,*)  (ZZZ(I),I=1,NATSOL)
      WRITE (6,578)(ZZZ(I),I=1,NATSOL)
      READ (51,*)
      READ (51,*) (DA(I),I=1,3)
c      WRITE (6,113) DA(1),DA(2),DA(3)
      READ (51,*)
      READ (51,*) SPC
c     write(*,*) 'SPC=',SPC

      IF(SPC.EQ.1)THEN

        alpha=1.00D0
        alpha2=0.00D0

        ALFA1 = alpha
        ALFA2 = alpha2
      ELSE

        alpha=0.7439762D0
        alpha2=0.1280119D0

        ALFA1 = alpha
        ALFA2 = alpha2
      ENDIF

      READ (51,*)
      READ (51,*) ICLSTR
      READ (51,*)
      IF(ICLSTR.EQ.1)THEN
      READ (51,*) BXLGTH
      READ (51,*)
      READ (51,*)
      ELSE
      READ (51,*)
      READ (51,*) 
      READ (51,*) RHO
      ENDIF
      READ (51,*)
      READ (51,*) XFAX,RBUF
      READ (51,*)
      READ (51,*) IEWLD
      WRITE (6,115) IEWLD
      READ (51,*)
      READ (51,*) KMAX
      READ (51,*)
      READ (51,*) IUNID
      READ (51,*)
      READ (51,*) ELFIELD,D,NQ
C      READ (51,*)
c      READ (51,*)IDIPCOR,PMAX,PZMAX 


      IF(ICLSTR.EQ.1)THEN
      WRITE(6,114)BXLGTH
      ELSE
      NWAT=NNAT(1+NSPECQ)
      BXLGTH=(DBLE(NWAT)/RHO)**THIRD3
      
      IF(ELFIELD.EQ.1) THEN
      BXLGTH = NQ*0.2D0*D
      ENDIF

      WRITE (6,114)BXLGTH
      ENDIF
      IF((ICLSTR.EQ.1).AND.(IEWLD.EQ.1))THEN
      WRITE(6,*)'Revisar constraint.dat,cluster con Ewald!'
      STOP
      ENDIF

55    FORMAT (/,2X,'ICON:',I3,2X,'NTIME:',I9,3X,'NSCAL:',I4,3X,
     &'NSAVE:',I9)
59    FORMAT (/,2X,'NDIP',I8,2X,'IRDF:',I8,3X,'ICMOT:',I4,/,3X,
     &'NEWV:',I9,'IPR:',I9,2X,'NANG:',I6 )
56    FORMAT (/,2X,'MASS OF THE PARTICLE:',G14.6)
57    FORMAT (/,2X,'TEMP:',G14.6,2X,'TEMPQ:',G10.5,'DELTAT:',G14.6)
60    FORMAT (/,2X,'MASA TERMOST:',G14.6,2X,'MASA TERMOQ:',G14.6)
58    FORMAT (/,2X,'SIGMA:',G14.6)
558   FORMAT (/,2X,'EPSILON:',G14.6)
578   FORMAT (/,2X,'ZZZ:',G14.6)
70    FORMAT (/,2X,'NUMERO DE ESPECIES CLASICAS:',I6)
656   FORMAT (/,2X,'NUMERO DE ATOMOS POR MOLEC.SOLVENTE: ',7I5)
111   FORMAT(/,2X,'SIGMA(1):',G14.6)
112   FORMAT(/,2X,'EPSILON(1):',G14.6)
113   FORMAT(/,2X,'DISTANCIAS INTRAMOLECULARES',G14.6,2X,G14.6,2X,G14.6)
114   FORMAT(/,2X,'TAMANIO DE LA CAJA (A)',G14.6)
115   FORMAT(/,2X,'SI IEWLD=1 LLAMA FURIM',I6)
116   FORMAT(/,2X,'RCUTOFF=',G14.6)
117   FORMAT(/,2X,'SISTEMA: CLUSTER')
118   FORMAT(/,2X,'SISTEMA: FASE CONDENSADA') 
119   FORMAT(/,2X,'DENSIDAD:',G14.6,2X)
120   FORMAT(/,2X,'AGUA PURA, SIN SOLUTO')      
1000  FORMAT(/,2X,'SOLUTO: ',20A5)
      IF(ICMOT.EQ.2.AND.IQMOT.NE.2)THEN
      WRITE(6,*)'OPTIMIZA SOLO LOS NUCLEOS'
      ELSEIF(IQMOT.EQ.2.AND.ICMOT.NE.2)THEN
      WRITE(6,*)'OPTIMIZA SOLO Q PARCIALES CLASICAS'
      PAUSE
      ENDIF
      
      IF(NSPECQ.EQ.0)WRITE(6,120)
      IF (ICLSTR.EQ.1)THEN
      WRITE(6,117)
      ELSE
      WRITE (6,118)
      WRITE (6,*)'BULK'
      ENDIF

      RETURN
      END




      
