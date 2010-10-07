      SUBROUTINE FURIM
      INCLUDE 'param'
        INCLUDE 'COMM'
      DIMENSION ELR(NAT,200),ELI(NAT,200),EMR(NAT,200)
      DIMENSION EMI(NAT,200),ENR(NAT,200),ENI(NAT,200)
      DIMENSION FURX(NAT),FURY(NAT),FURZ(NAT),FURQ(NAT)
      DIMENSION EXPIKR(NAT),EXPIKI(NAT)
      DIMENSION CLM(NAT),SLM(NAT)




      IF(3*NWAT.NE.NPART)THEN
      WRITE(*,*)'NO ESTA PREPARADO PARA EWALD CON SOLUTO'
      STOP
      ENDIF
      DO 10 I=1,NPART
       FURX(I)=ZERO
       FURY(I)=ZERO
       FURZ(I)=ZERO
       FURQ(I)=ZERO
10    CONTINUE

*-------------------Construye exp(iK.R) -------------------*
*------Calcula explicitamente para nx, ny, nz=0,1,   ------*

      ENEFUR=ZERO
      
      DO 20 I=1, NPART
      ELR(I,1)= ONE
      ELI(I,1)= ZERO
      EMR(I,1)= ONE
      EMI(I,1)= ZERO
      ENR(I,1)= ONE
      ENI(I,1)= ZERO
      ARGX=TWOPI*XX(I)
      ARGY=TWOPI*YY(I)
      ARGZ=TWOPI*ZZ(I)
      ELR(I,2)=COS(ARGX)
      EMR(I,2)=COS(ARGY)
      ENR(I,2)=COS(ARGZ)
      ELI(I,2)=SIN(ARGX)
      EMI(I,2)=SIN(ARGY)
      ENI(I,2)=SIN(ARGZ)
20    CONTINUE


*-------Los restantes exp(iK.R) por recurrencia ----------*
      
      DO 30 LM1=2,KMAX
      L=LM1+1
      
      DO 30 I=1,NPART
      ELR(I,L)=ELR(I,LM1)*ELR(I,2) - ELI(I,LM1)*ELI(I,2)
      ELI(I,L)=ELR(I,LM1)*ELI(I,2) + ELI(I,LM1)*ELR(I,2)
      EMR(I,L)=EMR(I,LM1)*EMR(I,2) - EMI(I,LM1)*EMI(I,2)
      EMI(I,L)=EMR(I,LM1)*EMI(I,2) + EMI(I,LM1)*EMR(I,2)
      ENR(I,L)=ENR(I,LM1)*ENR(I,2) - ENI(I,LM1)*ENI(I,2)
      ENI(I,L)=ENR(I,LM1)*ENI(I,2) + ENI(I,LM1)*ENR(I,2)
30    CONTINUE

*---------Suma sobre todos los vectores y particulas---------*
*--       Loop sobre particulas esta dentro del loop sobre K's       

      MKK = 0
      QPEFAC = ONE
      
      DO 40 L = 0, KMAX
      L1    = L + 1
      RL    = TWOPI*DBLE(L)
      RKL  = AXI * RL
      QFRFAC= -TWO * QPEFAC

      DO 50 M = -KMAX, KMAX
      RM   = TWOPI*DBLE(M)
      RKM =  RM*BYI
      RKSQ = RKL*RKL + RKM*RKM 
      IF (L**2+M**2.GT.KSQMAX) GO TO 50
      M1   = IABS(M)+1
      DSM  = DBLE(SIGN(1,M))
      
      DO 500 I = 1, NPART
      EMII   = DSM*EMI(I,M1)
      CLM(I) = EMR(I,M1)*ELR(I,L1) - EMII*ELI(I,L1)
      SLM(I) = EMR(I,M1)*ELI(I,L1) + EMII*ELR(I,L1)
500   CONTINUE
      
      DO 60 N = -KMAX, KMAX
      RN    = TWOPI*DBLE(N)
      RKN  =  RN*CZI
      RKSQ  = RKL*RKL + RKM*RKM + RKN*RKN
      IF(L**2+M**2+N**2.GT.KSQMAX.OR.RKSQ.LT.1.0D-12)GOTO 60
      RRKSQ = ONE/RKSQ
      N1    = IABS(N)+1
      DSN   = DBLE(SIGN(1,N))
      MKK   = MKK + 1
      AKN   = QA*EXP(-QB*RKSQ)*RRKSQ
      SUMR  = ZERO
      SUMI  = ZERO

      DO 80 I = 1, NPART
      ENII = DSN * ENI(I,N1)
c--- 
c      EXPIKR(I) = PC(I)* (CLM(I)*ENR(I,N1) - ENII*SLM(I))
c      EXPIKI(I) = PC(I)* (SLM(I)*ENR(I,N1) + ENII*CLM(I))
c      SUMR = SUMR + EXPIKR(I)
c      SUMI = SUMI + EXPIKI(I)
c--- 
      EXPIKR(I) =  (CLM(I)*ENR(I,N1) - ENII*SLM(I))
      EXPIKI(I) =  (SLM(I)*ENR(I,N1) + ENII*CLM(I))
      SUMR = SUMR + EXPIKR(I)*PC(I)
      SUMI = SUMI + EXPIKI(I)*PC(I)
80    CONTINUE
      

      QENEF = QPEFAC * AKN * (SUMR*SUMR + SUMI*SUMI)


C-----  ESTA ES LA ENERGIA  -----C

      ENEFUR = ENEFUR + QENEF
      QFACR  = QFRFAC * AKN * SUMR
      QFACI  = QFRFAC * AKN * SUMI

      DO 100 I=1,NPART

      QFORCE  = QFACI*EXPIKR(I)*PC(I) - EXPIKI(I)*QFACR*PC(I)
      QFORCEQ = QFACR*EXPIKR(I)      + QFACI*EXPIKI(I)

      FURX(I) = FURX(I) + RKL*QFORCE
      FURY(I) = FURY(I) + RKM*QFORCE
      FURZ(I) = FURZ(I) + RKN*QFORCE
      FURQ(I) = FURQ(I) + QFORCEQ

100   CONTINUE
      
60    CONTINUE
50    CONTINUE
      
         QPEFAC = TWO
   
40    CONTINUE   

c---- Si las cargas son fluctuantes, hay que agregar el termino
c---- 5.24 (pag.160 Allen & Tildesley) en la energia: EAUTO + EFF

c      if(one.ne.two)goto 777

      IF(IFLUC.EQ.1)THEN
      EAUTO = ZERO
      DO I=1,NPART
      EAUTO  =  EAUTO - PC(I)*PC(I)*RKAPPA/DSQRT(PI)
      FQ(I)  =  FQ(I) + TWO*PC(I)*RKAPPA/DSQRT(PI)
      ENDDO

      EFF = ZERO

c      if(one.ne.two)goto 771
c      App1=one/dsqrt(pi)- derf(da(1)/two)/da(1)
c      Bpp1=one/dsqrt(pi)- derf(da(3)/two)/da(3)
c      write(*,*)'A: ',app1
c      write(*,*)'B: ',bpp1
c      stop

      DO I = 1, NWAT
      AF1 = PC(I)  *  PC(I+NWAT)   * DERF(RKAPPA*DA(1))/DA(1)
      AF2 = PC(I)  *PC(I+2*NWAT)   * DERF(RKAPPA*DA(2))/DA(2)
      AF3 = PC(I+2*NWAT)*PC(I+NWAT)* DERF(RKAPPA*DA(3))/DA(3)

      EFF =  EFF - AF1 - AF2 - AF3

      FQ(I) = FQ(I) +               (AF1+AF2)/PC(I)
      FQ(I+NWAT) = FQ(I+NWAT) +     (AF1+AF3)/PC(I+NWAT)
      FQ(I+2*NWAT) = FQ(I+2*NWAT) + (AF2+AF3)/PC(I+2*NWAT)

      ENDDO

771   continue

      ENEFUR = ENEFUR + EAUTO + EFF

      ENDIF

777   continue

      NOFSET = 0
      DO 130 I = 1, NSPEC
      IF (I.GT.1) NOFSET = NOFSET+NNAT(I-1)
      DO 130 J = 1+NOFSET, NNAT(I)+NOFSET 


       FQ(J) = FQ(J) + FURQ(J)


       IF (I.EQ.1) THEN

       FX(J) = FX(J) + FURX(J)*ALFA1
       FY(J) = FY(J) + FURY(J)*ALFA1
       FZ(J) = FZ(J) + FURZ(J)*ALFA1
       J1 = J + NWAT
       FX(J1) = FX(J1) + FURX(J)*ALFA2
       FY(J1) = FY(J1) + FURY(J)*ALFA2
       FZ(J1) = FZ(J1) + FURZ(J)*ALFA2
       J1 = J1 + NWAT
       FX(J1) = FX(J1) + FURX(J)*ALFA2
       FY(J1) = FY(J1) + FURY(J)*ALFA2
       FZ(J1) = FZ(J1) + FURZ(J)*ALFA2

       ELSE

        FX(J) = FX(J) + FURX(J)
        FY(J) = FY(J) + FURY(J)
        FZ(J) = FZ(J) + FURZ(J)

       ENDIF

      write(*,*) 'zapato'
  
130   CONTINUE



      RETURN
      END



