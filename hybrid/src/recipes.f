!**********************************************************************
!    This file contains routines adapted from 'Numerical Recipes, 
!    The Art of Scientific Computing' by W.H. Press, S.A. Teukolsky, 
!    W.T. Veterling and B.P. Flannery, Cambridge U.P. 1987 and 1992.
!**********************************************************************
! The routines contained in this file are:
!     SUBROUTINE FOUR1
!     SUBROUTINE POLINT
!     SUBROUTINE SPLINE
!     SUBROUTINE SPLINT
!**********************************************************************


      SUBROUTINE FOUR1(DATA,NN,ISIGN)
!**********************************************************************
! Discrete Fourier transform. Modified and converted to double 
! precision from same routine in Numerical Recipes.
!**********************************************************************
! Input:
!   complex*16 DATA(NN) : Function to be Fourier transformed
!   integer    NN       : Number of points. Must be a power of 2
!   integer    ISIGN    : ISIG=+1/-1 => Direct/inverse transform
! Output:
!   complex*16 DATA(NN) : The direct Fourier transform (ISIG=+1), or
!                         NN times the inverse Fourier transf (ISIG=-1)
!**********************************************************************
      IMPLICIT NONE
      INTEGER          :: NN, ISIGN
      DOUBLE PRECISION :: DATA(2*NN)

      INTEGER          :: I, ISTEP, J, M, MMAX, N
      DOUBLE PRECISION :: TEMPI, TEMPR, THETA, WI, WPI, WPR, WR, WTEMP
      DOUBLE PRECISION, PARAMETER :: TWOPI=6.28318530717959D0,
     .  HALF=0.5D0, ONE=1.D0, TWO=2.D0, ZERO=0.D0

      N=2*NN
      J=1
      DO I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
        DO ! until following condition is met
          IF ((M.LT.2).OR.(J.LE.M)) EXIT
          J=J-M
          M=M/2
        END DO
        J=J+M
      END DO ! I
      MMAX=2
      DO ! until following condition is met
        IF (N.LE.MMAX) EXIT
        ISTEP=2*MMAX
        THETA=TWOPI/(ISIGN*MMAX)
        WPR=(-TWO)*SIN(HALF*THETA)**2
        WPI=SIN(THETA)
        WR=ONE
        WI=ZERO
        DO M=1,MMAX,2
          DO I=M,N,ISTEP
            J=I+MMAX
            TEMPR=WR*DATA(J)-WI*DATA(J+1)
            TEMPI=WR*DATA(J+1)+WI*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
          END DO ! I
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
        END DO ! M
        MMAX=ISTEP
      END DO ! until (N.LE.MMAX)

      END SUBROUTINE FOUR1



      SUBROUTINE POLINT(XA,YA,N,X,Y,DY) 
!*****************************************************************
! Polinomic interpolation. Modified and adapted to double 
! precision from same routine of Numerical Recipes.
! D. Sanchez-Portal, Oct. 1996
!*****************************************************************
! Input:
!   real*8  XA(N) : x values of the function y(x) to interpolate
!   real*8  YA(N) : y values of the function y(x) to interpolate
!   integer N     : Number of data points
!   real*8  X     : x value at which the interpolation is desired
! Output:
!   real*8  Y     : interpolated value of y(x) at X
!   real*8  DY    : accuracy estimate
!*****************************************************************

      IMPLICIT NONE
      INTEGER          :: N
      DOUBLE PRECISION :: XA(N),YA(N), X, Y, DY

      INTEGER          :: I, M, NS
      DOUBLE PRECISION :: C(N), D(N), DEN, DIF, DIFT, HO, HP, W
      DOUBLE PRECISION, PARAMETER :: ZERO=0.D0

      NS=1
      DIF=ABS(X-XA(1))
      DO I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
      END DO ! I
      Y=YA(NS)
      NS=NS-1
      DO M=1,N-1
        DO I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF (DEN.EQ.ZERO) STOP 'polint: ERROR. Two XAs are equal'
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
        END DO ! I
        IF (2*NS.LT.N-M) THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
      END DO ! M

      END SUBROUTINE POLINT



      SUBROUTINE SPLINE(DX,Y,N,YP1,YPN,Y2) 
!*********************************************************** 
! Cubic Spline Interpolation. Adapted from Numerical Recipes 
! routine of same name for a uniform grid and double precision
! D. Sanchez-Portal, Oct. 1996.
! Input:
!   real*8  DX   : x interval between data points
!   real*8  Y(N) : value of y(x) at data points
!   integer N    : number of data points
!   real*8  YP1  : value of dy/dx at X1 (first point)
!   real*8  YPN  : value of dy/dx at XN (last point)
! Output:
!   real*8  Y2(N): array to be used by routine SPLINT
! Behavior:
! - If YP1 or YPN are larger than 1E30, the natural spline
!   condition (d2y/dx2=0) at the corresponding edge point.
!************************************************************

      IMPLICIT NONE
      INTEGER          :: N
      DOUBLE PRECISION :: DX, Y(N), YP1, YPN, Y2(N)

      INTEGER          :: I, K
      DOUBLE PRECISION :: QN, P, SIG, U(N), UN
      DOUBLE PRECISION, PARAMETER :: YPMAX=0.99D30, 
     .  HALF=0.5D0, ONE=1.D0, THREE=3.D0, TWO=2.D0, ZERO=0.D0
    
      IF (YP1.GT.YPMAX) THEN
        Y2(1)=ZERO
        U(1)=ZERO
      ELSE
        Y2(1)=-HALF
        U(1)=(THREE/DX)*((Y(2)-Y(1))/DX-YP1)
      ENDIF
      DO I=2,N-1
        SIG=HALF
        P=SIG*Y2(I-1)+TWO
        Y2(I)=(SIG-ONE)/P
        U(I)=(THREE*( Y(I+1)+Y(I-1)-TWO*Y(I) )/(DX*DX)
     .       -SIG*U(I-1))/P
      END DO ! I
      IF (YPN.GT.YPMAX) THEN
        QN=ZERO
        UN=ZERO
      ELSE
        QN=HALF
        UN=(THREE/DX)*(YPN-(Y(N)-Y(N-1))/DX)
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+ONE)
      DO K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
      END DO ! K

      END SUBROUTINE SPLINE



      SUBROUTINE SPLINT(DX,YA,Y2A,N,X,Y,DYDX) 
!***************************************************************
! Cubic Spline Interpolation. Adapted from Numerical Recipes 
! routine of same name for a uniform grid, double precision,
! and to return the function derivative in addition to its value
! D. Sanchez-Portal, Oct. 1996.
! Input:
!   real*8  DX    : x interval between data points
!   real*8  YA(N) : value of y(x) at data points
!   real*8  Y2A(N): array returned by routine SPLINE
!   integer N     : number of data points
!   real*8  X     : point at which interpolation is desired
!   real*8  Y     : interpolated value of y(x) at point X
!   real*8  DYDX  : interpolated value of dy/dx at point X
!***************************************************************

      IMPLICIT NONE
      INTEGER          :: N
      DOUBLE PRECISION :: DX, YA(N), Y2A(N), X, Y, DYDX

      INTEGER          :: NHI, NLO
      DOUBLE PRECISION :: A, B
      DOUBLE PRECISION, PARAMETER ::
     .    ONE=1.D0, THREE=3.D0, SIX=6.D0, ZERO=0.D0

      IF (DX.EQ.ZERO) STOP 'splint: ERROR: DX=0'
      NLO=INT(X/DX)+1
      NHI=NLO+1
      A=NHI-X/DX-1
      B=ONE-A
      Y=A*YA(NLO)+B*YA(NHI)+
     .  ((A**3-A)*Y2A(NLO)+(B**3-B)*Y2A(NHI))*(DX**2)/SIX
      DYDX=(YA(NHI)-YA(NLO))/DX+
     .     (-((THREE*(A**2)-ONE)*Y2A(NLO))+
     .     (THREE*(B**2)-ONE)*Y2A(NHI))*DX/SIX

      END SUBROUTINE SPLINT
