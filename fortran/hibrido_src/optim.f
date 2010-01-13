      SUBROUTINE OPTIM(NATSOL)
      INCLUDE 'COMM'
      DIMENSION AXX(NAT),AYY(NAT),AZZ(NAT)
 

      AA=ZERO
c     AA=0.001D0
      NOFSET = 0
      DO 100 I = 1, NSPECQ+NATSOL
      IF (I.GT.1) NOFSET = NOFSET + NNAT(I-1)
      DO 100 J= 1+NOFSET,NNAT(I)+NOFSET 
      XX(J) = X(J) + (X(J)-X0(J))*AA + FFF(I) * FX(J)
      YY(J) = Y(J) + (Y(J)-Y0(J))*AA + FFF(I) * FY(J)
      ZZ(J) = Z(J) + (Z(J)-Z0(J))*AA + FFF(I) * FZ(J)
 
c      write(*,*)'x OPT ',j,fx(j),fx(j)*fff(i),x(j),xx(j)
c      write(*,*)'y OPT ',j,fy(j),fy(j)*fff(i),y(j),yy(j)
c      write(*,*)'z OPT ',j,fz(j),fz(j)*fff(i),z(j),zz(j)
     
c      AXX(J) = XX(J)
c      AYY(J) = YY(J)
c      AZZ(J) = ZZ(J)
100   CONTINUE


      CALL GAMMA(NATSOL) 
c      do i=1,npart
c      write(*,*)'opt-gam ',i,xx(i),yy(i),zz(i)
c      enddo



     
C-----  UPDATE POSITIONS

      DO 12 I = 1, NPART

      AX1 = X0(I)
      X0(I) = X(I)
      X(I) = XX(I)
      XX(I) = AX1

      AY1= Y0(I)
      Y0(I) = Y(I)
      Y(I) = YY(I)
      YY(I) = AY1

      AZ1= Z0(I)
      Z0(I) = Z(I)
      Z(I) = ZZ(I)
      ZZ(I) = AZ1
c      write(*,*)'opt ',i,x(i),y(i),z(i)
 12   CONTINUE


      RETURN
      END


