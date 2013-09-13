      SUBROUTINE OPTIMQ(NATSOL)
      use hibrido_common
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      AA=ZERO
c      AA=0.001     
      DO 80 I=1+NATOM,NWAT+NATOM
      I1=I+NWAT
      I2=I1+NWAT
      PCC(I) =PC(I)+(PC(I)-PC0(I))*AA + DEL/WWQ*FQ(I)
      PCC(I1) =PC(I1)+(PC(I1)-PC0(I2))*AA + DEL/WWQ*FQ(I1)
      PCC(I2) =PC(I2)+(PC(I2)-PC0(I1))*AA + DEL/WWQ*FQ(I2)

*----Constraint c/molecula neutra
      SUMQ=THIRD3*(PCC(I)+PCC(I1)+PCC(I2))
      PCC(I) = PCC(I) -SUMQ
      PCC(I1)= PCC(I1)-SUMQ
      PCC(I2)= PCC(I2)-SUMQ

      SUMQ=PCC(I)+PCC(I1)+PCC(I2)
      IF(DABS(SUMQ/EE).GT.TOLL)THEN
      WRITE(6,*)'PROBLEMAS CONSTRAINTS QS (optimq.f), EN PASO: ',ITEL
      WRITE(6,*)'I,SUMQ/E,CARGAS',I,SUMQ/EE,PCC(I)/EE,
     &  PCC(I1)/EE,PCC(I2)/EE
      STOP
      ENDIF






80    CONTINUE
      
 

C-----UPDATE CHARGES 

      DO 12 I = 1+NATOM, NPART
      AQ1 = PC0(I)
      PC0(I) = PC(I)
      PC(I) = PCC(I)
      PCC(I) = AQ1
12    CONTINUE

      RETURN
      END


