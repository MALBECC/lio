          SUBROUTINE PRINTER (RT,HISTO,HISTO2,DELP,DELPZ)
	  use hibrido_common
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          DIMENSION HISTO(100,100),HISTO2(100)

          SUMA1 = ZERO
          SUMA2 = ZERO
          SUMA3 = ZERO
          SUMA4 = ZERO
          SUMA5 = ZERO
          SUMA6 = ZERO
          SUMA7 = ZERO
          SUMA8 = ZERO

	  DO 1 I = 1, IGGRID

	  R = DELR*DBLE(I)
	  RMIN = R - DELR/TWO
	  RMAX = R + DELR/TWO
	  IF (ICLSTR.EQ.1) THEN
	  VOL=FPI*(RMAX**3-RMIN**3)/THREE*DBLE(NWAT)
	  ELSE
	  VOL=FPI*(RMAX**3-RMIN**3)/THREE*RHO*DBLE(NWAT)
	  ENDIF
           
          GR1 =  DBLE(KG(1,I))/VOL*RT
          GR2 =  DBLE(KG(2,I))/VOL/TWO*RT
          GR3 =  DBLE(KG(3,I))/VOL/FOUR*RT
          GR4 =  DBLE(KG(4,I))/VOL*RT
          GR5 =  DBLE(KG(5,I))/VOL*RT
          GR6 =  DBLE(KG(6,I))/VOL/TWO*RT
          GR7 =  DBLE(KG(7,I))/VOL/TWO*RT
          GR8 =  DBLE(KG(8,I))/VOL/TWO*RT

          SUMA1 = SUMA1 + DBLE(KG(1,I))/DBLE(NWAT)*RT
          SUMA2 = SUMA2 + DBLE(KG(2,I))/DBLE(NWAT)*RT
          SUMA3 = SUMA3 + DBLE(KG(3,I))/DBLE(2*NWAT)*RT
          SUMA4 = SUMA4 + DBLE(KG(4,I))/DBLE(NWAT)*RT
          SUMA5 = SUMA5 + DBLE(KG(5,I))/DBLE(NWAT)*RT
          SUMA6 = SUMA6 + DBLE(KG(6,I))/DBLE(2*NWAT)*RT
          SUMA7 = SUMA7 + DBLE(KG(7,I))/DBLE(2*NWAT)*RT
          SUMA8 = SUMA8 + DBLE(KG(8,I))/DBLE(2*NWAT)*RT
         
          WRITE (16,552)R,GR1,GR2,GR3,SUMA1,SUMA2,SUMA3
          WRITE (17,552)R,GR4,GR5,GR6,GR7,GR8,SUMA4,SUMA5,SUMA6,SUMA7,
     &    SUMA8
552       FORMAT(2X,11G16.7)

1         CONTINUE 


C-----IMPRIME EN FILE 11: fort.11
       SUMQK  = ZERO
       SUMGR1 = ZERO
       SUMGR2 = ZERO

       DO 11 I = 1, IGRIL+1
       DO 12 J = 1, IGRIL+1
      
       Z3 = DELR2*(DBLE(I)-GRIL/TWO-ONE)
       R3 = DELR2*(DBLE(J)-ONE)

       IF(R3.EQ.ZERO)THEN
       QKK = QK(I,J)*RT
       GRR1= DBLE(KG1(1,I,J))*RT
       GRR2= DBLE(KG1(2,I,J))/FOUR*RT
       GOTO 45

       ENDIF

       RMIN = R3 - DELR2/TWO
       RMAX = R3 + DELR2/TWO
       VOL = PI*(RMAX*RMAX-RMIN*RMIN)*DBLE(NWAT)

       QKK = QK(I,J)/VOL*DBLE(NWAT)*RT
       GRR1= DBLE(KG1(1,I,J))/VOL*RT
       GRR2= DBLE(KG1(2,I,J))/VOL/FOUR*RT
       
45     CONTINUE
       SUMQK = SUMQK + QK(I,J)*RT
       SUMGR1 = SUMGR1 + DBLE(KG1(1,I,J))/DBLE(NWAT)*RT
       SUMGR2 = SUMGR2 + DBLE(KG1(2,I,J))/DBLE(2*NWAT)*RT
     
       WRITE(11,552)Z3,R3,QKK,GRR1,GRR2,SUMQK,SUMGR1,SUMGR2
12     CONTINUE
       WRITE(11,*)
11     CONTINUE

C---CORRELACION DIPOLAR
C---AL FINAL

c       DO II = 1, 100 
       DO J = 1, 100

          R = DBLE(J-1)*DELP
          RMIN = R - DELP/2.D0
          RMAX = R - DELP/2.D0
c          ZP = (DBLE(II)-50.D0-1.D0)*DELPZ
       VOL = FPI*(RMAX**3 - RMIN**3)/THREE*DBLE(NWAT)
       IF(VOL.EQ.0.D0)VOL=1.D0
c          HISTO(II,J) = HISTO(II,J)/VOL*RT
          HISTO2(J) = HISTO2(J)/VOL*RT
c       WRITE(12,54)R,ZP,HISTO(II,J)
       WRITE(12,54)R,HISTO2(J)
       ENDDO
c       WRITE(12,*)
c       ENDDO

54     FORMAT(1X,3G11.5)


      RETURN
      END


	  
