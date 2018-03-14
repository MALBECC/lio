c--------------------------------------------------------------
c this small program is taken from old version of gdfmol
c it is the implementation of the Obara-Saika method for
c the evaluation of F(m,T), using a 2 branch calculation
c DEBUGGING VERSION, this is the attempt to generalize
c and improve previous version ( up to F16 )
c Ref: JCP 84 3963 (1986)
c it seems to work
c This is the version that should be included in definitive
c program
c 11 March 1992 -----------------------------------------------
c
      FUNCTION FUNCT(N,T)                                               
      USE garcha_mod
      IMPLICIT REAL*8 (A-H,O-Z)                                         
!      COMMON /TABLE/ STR(880,0:21),FAC(0:16)
C
      if (T.lt.0.0D0) then
       write(*,*) 'Problems',T
       T=abs(T)
      endif
C                                                                       
      IF (T.LE.43.975D0)   THEN                          
       IT = 20.0D0 * (T + 0.025D0)                                        
       TI = DFLOAT(IT)                                                   
       IT=IT + 1                                                         
       DELT = T - 0.05D0 * TI                                             
       DELT3 = DELT * 0.333333333333333D0                                
       DELT4 = 0.25D0 * DELT                                             
       DELT2 = DELT4 + DELT4                                             
       DELT5 = 0.20D0 * DELT                                             
c
c
c
       TF0 = STR(IT,N)                                                    
       TF1 = STR(IT,N+1)                                                    
       TF2 = STR(IT,N+2)                                                    
       TF3 = STR(IT,N+3)                                                    
       TF4 = STR(IT,N+4)                                                    
       TF5 = STR(IT,N+5)                                                    
c
       FCAP=TF0-DELT*( TF1-DELT2*(TF2-DELT3*(TF3-DELT4*               
     >    (TF4-DELT5*TF5))))                                           
       FUNCT = FCAP                                                   
       RETURN                                                            
C                                                                       
       
       ELSE
      FUNCT=FAC(N)*1.D0/(T**N*dsqrt(T))
     
      ENDIF
c
      END            
  

C ----------------------------------------------                        
c subroutine for generating tables, used later on in the Taylor expansion
c for the incomplete Gamma functions
c-----------------------------------------------
      SUBROUTINE GENERF
      USE garcha_mod
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (sqpi=1.77245385090551588D0)
!      COMMON /TABLE/ STR(880,0:21),FAC(0:16)
c
c loop over T values in the table ( 0. to 43.95 , interval 0.05)
c
      DO 799 I=1,880
       T = 0.05D0 * DFLOAT(I-1)
       Y = DEXP(-T)
       U = T + T
       F = FMCH(22,T)
c
c loop over order of incomple Gamma functions ( 0 to 21, the ones
c necessary for evaluating orders 0-16)
c
      DO 99 M=21,0,-1
       W=2.D0*DFLOAT(M)+1.D0
       STR(I,M) = (Y + U * F)/W
       F=STR(I,M)
  99  CONTINUE
c
  799 CONTINUE
c calculation of the function [(v+1/2)/2**v+1*sqrt(pi) ]
c
      FAC(0)=sqpi/2.D0
      do n=1,16
       FAC(n)=FAC(n-1)*(2*n-1)/2
      enddo
c
      RETURN
      END

      SUBROUTINE GENERFS
      IMPLICIT REAL*4 (A-H,O-Z)
      parameter (sqpi=1.772453850905E0)
      COMMON /TABLES/ STRS(880,0:21),FACS(0:16)
c
c loop over T values in the table ( 0. to 43.95 , interval 0.05)
c
      DO 799 I=1,880
       T = 0.05E0 * FLOAT(I-1)
       Y = EXP(-T)
       U = T + T
       F = FMCHS(22,T)
c
c loop over order of incomple Gamma functions ( 0 to 21, the ones
c necessary for evaluating orders 0-16)
c
      DO 99 M=21,0,-1
       W=2.E0*FLOAT(M)+1.E0
       STRS(I,M) = (Y + U * F)/W
       F=STRS(I,M)
  99  CONTINUE
c
  799 CONTINUE
c calculation of the function [(v+1/2)/2**v+1*sqrt(pi) ]
c
      FACS(0)=sqpi/2.E0
      do n=1,16
       FACS(n)=FACS(n-1)*(2*n-1)/2
      enddo
c
      RETURN
      END

c-------------------
c same version as in old version of GDFMOL
      FUNCTION FMCH(M,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (sqpi=1.77245385090551588D0)
c
      Y=DEXP(-X)
      IF (X.GT.10.D0) GO TO 20
      A = DFLOAT(M)
      A = A + 0.5D0
      TERM = 1.0D0 / A
      PTLSUM = TERM
      DO 11 I=2,50
      TERM = TERM * X / (A + DFLOAT(I-1))
      PTLSUM = PTLSUM + TERM
      IF(TERM/PTLSUM.LT.1.0D-12)GO TO 12
   11 CONTINUE 
      STOP     
   12 FMCH = 0.5D0 * PTLSUM * Y
      RETURN   
c
   20 A = DFLOAT(M)
      B = A + 0.5D0
      XD = 1.0D0 / X
      APPROX = 0.886226925452758D0 * DSQRT(XD) * XD**M
      DO 22 I=1,M
   22 APPROX = APPROX * (B - DFLOAT(I))
      FIMULT = 0.5D0 * Y * XD
      PTLSUM = 0.0D0
      IF(FIMULT.EQ.0.0D0)GO TO 25
      FIPROP = FIMULT/APPROX
      TERM = 1.0D0
      PTLSUM = TERM
      NOTRMS = X
      NOTRMS = NOTRMS + M
      DO 24 I=2,NOTRMS
      A = B - DFLOAT(I-1)
      TERM = TERM * A * XD
      PTLSUM = PTLSUM + TERM
      IF(DABS(TERM*FIPROP/PTLSUM).LE.1.0D-10)GO TO 25
   24 CONTINUE
      WRITE (IW,999) M,X
   25 FMCH = APPROX - FIMULT * PTLSUM
      RETURN
  999 FORMAT ('23H0FMCH DOES NOT CONVERGE',I6,1PE16.8)
      END
C     -----------------------
      FUNCTION FMCHS(M,X)
      IMPLICIT REAL*4 (A-H,O-Z)
      parameter (sqpi=1.77245385090E0)
c
      Y=EXP(-X)
      IF (X.GT.10.E0) GO TO 20
      A = FLOAT(M)
      A = A + 0.5E0
      TERM = 1.0E0 / A
      PTLSUM = TERM
      DO 11 I=2,50
      TERM = TERM * X / (A + FLOAT(I-1))
      PTLSUM = PTLSUM + TERM
      IF(TERM/PTLSUM.LT.1.0E-12)GO TO 12
   11 CONTINUE
      STOP
   12 FMCHS = 0.5E0 * PTLSUM * Y
      RETURN
c
   20 A = FLOAT(M)
      B = A + 0.5E0
      XD = 1.0E0 / X
      APPROX = 0.886226925452758E0 * SQRT(XD) * XD**M
      DO 22 I=1,M
   22 APPROX = APPROX * (B - FLOAT(I))
      FIMULT = 0.5E0 * Y * XD
      PTLSUM = 0.0E0
      IF(FIMULT.EQ.0.0E0)GO TO 25
      FIPROP = FIMULT/APPROX
      TERM = 1.0E0
      PTLSUM = TERM
      NOTRMS = X
      NOTRMS = NOTRMS + M
      DO 24 I=2,NOTRMS
      A = B - FLOAT(I-1)
      TERM = TERM * A * XD
      PTLSUM = PTLSUM + TERM
      IF(ABS(TERM*FIPROP/PTLSUM).LE.1.0E-10)GO TO 25
   24 CONTINUE
      WRITE (IW,999) M,X
   25 FMCHS = APPROX - FIMULT * PTLSUM
      RETURN
  999 FORMAT (24H0FMCHS DOES NOT CONVERGE,I6,1PE16.8)
      END
C     -----------------------

