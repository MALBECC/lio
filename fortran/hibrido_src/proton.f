      SUBROUTINE PROTON(FXH1,FYH1,FZH1,FXH11,FYH11,FZH11,
     > FXH2,FYH2,FZH2,ff,FZ5,FZ1,FZ2)
      INCLUDE 'param'
        INCLUDE 'COMM'

      DIMENSION ff(NAT,3)
c------------
      IPR2=5 
c------------

c--Fuerza neta sobre H (1)
      FXH= FX(1)           
      FYH= FY(1)         
      FZH= FZ(1)         
c--Contribucion a F(H) proveniente de Eq-c:
c   1 ) Electrostatica entre Z y Q parciales (FXH1,FYH1,FZH1)
c   11) De la integral de Rho * Q parciales

      FXH11= -FXH11*HH/A0          
      FYH11= -FYH11*HH/A0         
      FZH11= -FZH11*HH/A0             
  
c--Contribucion a F(H) proveniente de Eq
      FXH2 = -FXH2 *HH/A0          
      FYH2 = -FYH2 *HH/A0         
      FZH2 = -FZH2 *HH/A0             

c--- Control:
      IF(NDFT.EQ.1)THEN
      P1 = FXH1 + FXH11 + FXH2
      P2 = FYH1 + FYH11 + FYH2
      P3 = FZH1 + FZH11 + FZH2
 
      IF(FXH-P1.GT.1.D-07)WRITE(*,*)'FXH.NE.P1',FXH,P1
      IF(FYH-P2.GT.1.D-07)WRITE(*,*)'FYH.NE.P2',FYH,P2
      IF(FZH-P3.GT.1.D-07)WRITE(*,*)'FZH.NE.P3',FZH,P3
      ENDIF

c--- EJE OXIG-OXIG= EJE Z5(de H2O a NO3)
c--- IA= OXIGENO(NO3); IB= OXIGENO(H2O)
      IA=4
      IB=6
      XE=X(IA)-X(IB)
      YE=Y(IA)-Y(IB)
      ZE=Z(IA)-Z(IB)
      EJEZ=DSQRT(XE**2+YE**2+ZE**2)
      XE=XE/EJEZ
      YE=YE/EJEZ
      ZE=ZE/EJEZ
     
C---  PROYECCION de FXH SOBRE EL EJE Z5
      FZ5 = FXH *XE + FYH *YE + FZH *ZE
      FZ11= FXH11*XE + FYH11*YE + FZH11*ZE
      FZ1 = FXH1*XE + FYH1*YE + FZH1*ZE
      FZ2 = FXH2*XE + FYH2*YE + FZH2*ZE
      
      WRITE(90,*)'IA=4  ,  IB=6'   
      IF(MOD(ITEL,IPR2).EQ.0)THEN
      WRITE(90,110)ITEL,FZ5,FZ11,FZ1,FZ2
      ENDIF
110   FORMAT(I10,2X,4G15.7)


      
      RETURN
      END
