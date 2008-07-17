      SUBROUTINE QTMOTB
      INCLUDE 'COMM'
      
      
      DO 80 I=NATOM+1,NATOM+NWAT
      I1=I+NWAT
      I2=I1+NWAT
      PCC(I)   = TWO*PC(I) -PC0(I) +DEL/WWQ*FQ(I)
      PCC(I1)  = TWO*PC(I1)-PC0(I1)+DEL/WWQ*FQ(I1)
      PCC(I2)  = TWO*PC(I2)-PC0(I2)+DEL/WWQ*FQ(I2)

      SUMQ=THIRD3*(PCC(I)+PCC(I1)+PCC(I2))
      PCC(I) = PCC(I) -SUMQ
      PCC(I1)= PCC(I1)-SUMQ
      PCC(I2)= PCC(I2)-SUMQ
      
      SUMQ=PCC(I)+PCC(I1)+PCC(I2)
      IF(DABS(SUMQ/EE).GT.TOLL)THEN
      WRITE(6,*)'PROBLEMAS CONSTRAINTS QS (qemot.f), EN PASO: ',ITEL
      WRITE(6,*)'I,SUMQ/E,CARGAS',I,SUMQ/EE,PCC(I)/EE,
     &  PCC(I1)/EE,PCC(I2)/EE
      STOP
      ENDIF


80    CONTINUE
 
 
      DO I=NATOM+1,NATOM+NWAT
      I1=I+NWAT
      I2=I1+NWAT

      VQ(I) =FACTV * (PCC(I)-PC0(I))
      VQ(I1)=FACTV * (PCC(I1)-PC0(I1))
      VQ(I2)=FACTV * (PCC(I2)-PC0(I2))

      VQCM=VQ(I)+VQ(I1)+VQ(I2)    
      VQCM=VQCM*THIRD3

      VQ(I)=VQ(I)-VQCM
      VQ(I1)=VQ(I1)-VQCM
      VQ(I2)=VQ(I2)-VQCM

      VQ0(I) =VQ(I)
      VQ0(I1)=VQ(I1)
      VQ0(I2)=VQ(I2)

      ENDDO
   
      TEMPOL = ZERO
      
      DO I = NATOM+1, NPART
     
      TEMPOL=TEMPOL+VQ0(I)*VQ0(I)*WWQ

      ENDDO
c lau agregado
C=========temperatura constante
C=========Solvente-cargas
      tesq= 1. + (DELTAT/QSQ)*((TEMPRQQ/(TEMPOL/GDQ)) -1.)
c      write (*,*) 'lau', TEMPOL/GDQ,TEMPRQQ,tesq
      tesq=sqrt(tesq)

C-----  UPDATE CHARGES 

      DO 12 I = NATOM+1, NPART
      AQ1 = PC0(I)
      PC0(I) = PC(I)
      PC(I) =tesq*(PCC(I)-PC(I))+PC(I)
      PCC(I) = AQ1

c  lau agregado


12    CONTINUE

      RETURN
      END
