!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!  This small program is taken from old version of gdfmol
!  it is the implementation of the Obara-Saika method for
!  the evaluation of F(m,T), using a 2 branch calculation
!  DEBUGGING VERSION, this is the attempt to generalize
!  and improve previous version ( up to F16 ).
!  Ref: JCP 84 3963 (1986)
!  it seems to work
!
!  This is the version that should be included in definitive
!  program
!  11 March 1992
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      FUNCTION FUNCT(N,T)
      USE garcha_mod
      IMPLICIT REAL*8 (A-H,O-Z)
      integer, intent(in)    :: N
      real*8 , intent(inout) :: T
      integer                :: IT
!
      if (T.lt.0.0D0) then
       write(*,*) 'Problems',T
       T=abs(T)
      endif
!
      IF (T.LE.43.975D0) THEN
       IT = 20.0D0 * (T + 0.025D0)
       TI = DFLOAT(IT)
       IT=IT + 1
       DELT = T - 0.05D0 * TI
       DELT3 = DELT * 0.333333333333333D0
       DELT4 = 0.25D0 * DELT
       DELT2 = DELT4 + DELT4
       DELT5 = 0.20D0 * DELT
!
       TF0 = STR(IT,N)
       TF1 = STR(IT,N+1)
       TF2 = STR(IT,N+2)
       TF3 = STR(IT,N+3)
       TF4 = STR(IT,N+4)
       TF5 = STR(IT,N+5)
!
       FCAP=TF0-DELT*(TF1-DELT2*(TF2-DELT3*(TF3-DELT4*(TF4-DELT5*TF5))))
       FUNCT = FCAP
       RETURN
!
      ELSE
         FUNCT=FAC(N)*1.D0/(T**N*dsqrt(T))

      ENDIF
!
      END FUNCTION FUNCT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
