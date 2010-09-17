      SUBROUTINE TEPQ
      INCLUDE 'param'
        INCLUDE 'COMM'
      
*-----Velocid. al azar de 1 distr. gausseana
      TEMPOL = ZERO
      DO  1 J=NATOM+1,NATOM+NWAT
      J1=J+NWAT
      J2=J+2*NWAT
      VQ(J) =VFQ*GAUSSN()
      VQ(J1)=VFQ*GAUSSN()
      VQ(J2)=VFQ*GAUSSN()
      TEMPOL = TEMPOL + WWQ*(VQ(J)**2+VQ(J1)**2+VQ(J2)**2)
1     CONTINUE
      TEMPOL1 = TEMPOL/GDQ

*-----Le resto VCM a todas las cargas
*-----(en cada molecula) 
      VQCM= ZERO
      DO I=NATOM+01,NATOM+NWAT
      I1=I+NWAT
      I2=I+2*NWAT
      VQCM=(VQ(I)+VQ(I1)+VQ(I2))*THIRD
      VQ(I) =VQ(I) -VQCM
      VQ(I1)=VQ(I1)-VQCM
      VQ(I2)=VQ(I2)-VQCM
      ENDDO


      DO 20 I=NATOM+1,NATOM+NWAT
      I1=I +NWAT
      I2=I1+NWAT

*-----Un paso atras
      PC(I) = PC0(I)
      PC(I1) = PC0(I1)
      PC(I2) = PC0(I2)
      VQ0(I)=VQ(I)
      VQ0(I1)=VQ(I1)
      VQ0(I2)=VQ(I2)

*-----Cargas. nuevas
      PCC(I)=PC(I)+VQ0(I)*DELTAT
      PCC(I1)=PC(I1)+VQ0(I1)*DELTAT
      PCC(I2)=PC(I2)+VQ0(I2)*DELTAT

*-----Constraint c/molecula neutral
      SUMQ=THIRD3*(PCC(I)+PCC(I1)+PCC(I2))
      PCC(I) = PCC(I) -SUMQ
      PCC(I1)= PCC(I1)-SUMQ
      PCC(I2)= PCC(I2)-SUMQ

      SUMQ=PCC(I)+PCC(I1)+PCC(I2)
      IF(DABS(SUMQ/EE).GT.TOLL)THEN
      WRITE(6,*)'PROBLEMAS CONSTRAINTS QS (tepq.f), EN PASO: ',ITEL
      WRITE(6,*)'I,SUMQ/E,CARGAS',I,SUMQ/EE,PCC(I)/EE,
     &  PCC(I1)/EE,PCC(I2)/EE
      STOP
      ENDIF

20    CONTINUE
      

*-----Temperatura media

      TEMPOL=ZERO

      DO 35 J=NATOM+1,NPART
      PC(J)=PCC(J)
      VQ0(J) = (PC(J)-PC0(J))/DELTAT
      VQ(J) = VQ0(J)
      VQ2=VQ0(J)*VQ0(J)
      TEMPOL=TEMPOL+VQ2*WWQ
35    CONTINUE

      TEMPPQ=TEMPOL/GDQ
c      write(*,*)'tepq ',tempol,gdq

*--------TEMPAV1: despues de randomizar--------------*
*        TEMPP  : despues de llamar gamma y          *
*		  recalcular posiciones y velocids---*

      WRITE(6,19) TEMPRQQ,TEMPOL1
19     FORMAT (/,2X, 'TERMOSTATO:, TEMP1 ',G14.6,2X,
     & G14.6,2X,G14.6)
      WRITE (6,*)itel, ' ENTRE A TEPQ y TEMP ES:', TEMPPQ

      RETURN
      END

