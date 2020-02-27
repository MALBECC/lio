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
subroutine init_math
   use math_data, only: STR, FAC
   implicit none
   integer :: icount, jcount
   double precision :: T, Y, U, F, W
   double precision, parameter :: SQPI = 1.77245385090551588D0

   ! Loop over T values in the table ( 0. to 43.95 , interval 0.05)
   do icount = 1, 880
      T = 0.05D0 * DFLOAT(icount -1)
      Y = DEXP(-T)
      U = T + T
      F = FMCH(22,T)

      ! loop over order of incomple Gamma functions ( 0 to 21, the ones
      ! necessary for evaluating orders 0-16)
      do jcount = 21, 0, -1
         W = 2.0D0 * DFLOAT(jcount) + 1.0D0
         STR(icount,jcount) = (Y + U * F) / W
         F = STR(icount,jcount)
      enddo
   enddo

   ! Calculation of the function [(v+1/2)/2**v+1*sqrt(pi) ]
   FAC(0) = SQPI / 2.0D0
   do icount = 1, 16
      FAC(icount) = FAC(icount -1) * (2 * icount -1) / 2
   enddo

end subroutine init_math

function FMCH(M,X)
   implicit none
   integer         , intent(in) :: M
   double precision, intent(in) :: X

   double precision :: Y, A, PTLSUM, TERM, B, XD, APPROX, FIMULT, FMCH, FIPROP,&
                       NOTRMS
   integer          :: I

   Y = DEXP(-X)

   if (.not. (X .gt. 10.0D0)) then
      A = DFLOAT(M) + 0.5D0
      TERM   = 1.0D0 / A
      PTLSUM = TERM
      do I = 2, 50
         TERM = TERM * X / (A + DFLOAT(I-1))
         PTLSUM = PTLSUM + TERM
         if ((TERM / PTLSUM) .LT. 1.0D-12) exit
      enddo
      FMCH = 0.5D0 * PTLSUM * Y
      return
   else
      A      = DFLOAT(M)
      B      = A + 0.5D0
      XD     = 1.0D0 / X
      APPROX = 0.886226925452758D0 * DSQRT(XD) * XD**M

      do I = 1, M
         APPROX = APPROX * (B - DFLOAT(I))
      enddo

      FIMULT = 0.5D0 * Y * XD
      PTLSUM = 0.0D0
      if (FIMULT .eq. 0.0D0 ) then
         FMCH = APPROX
         return
      endif

      FIPROP = FIMULT / APPROX
      TERM   = 1.0D0
      PTLSUM = TERM
      NOTRMS = X
      NOTRMS = NOTRMS + M

      do I = 2, int(NOTRMS)
         A      = B - DFLOAT(I-1)
         TERM   = TERM * A * XD
         PTLSUM = PTLSUM + TERM
         if (DABS(TERM * FIPROP / PTLSUM) .le. 1.0D-10) exit
      enddo
      FMCH = APPROX - FIMULT * PTLSUM
   endif
return
end function FMCH

function funct(N, T)
   use math_data, only: STR, FAC
   implicit none
   integer         , intent(in)    :: N
   double precision, intent(inout) :: T

   double precision :: TI, DELT, DELT2, DELT3, DELT4, DELT5, TF0, TF1, TF2, &
                       TF3, TF4, TF5, funct, IT

   if (T .lt. 0.0D0) then
      write(*,'(A)') 'Problems with Boys function'
      T = abs(T)
   endif

   if (T .le. 43.975D0) then
      IT = 20.0D0 * (T + 0.025D0)
      TI = IT
      IT = IT + 1
      DELT = T - 0.05D0 * TI
      DELT3 = DELT * 0.333333333333333D0
      DELT4 = 0.25D0 * DELT
      DELT2 = DELT4 + DELT4
      DELT5 = 0.20D0 * DELT

      TF0 = STR(int(IT),N)
      TF1 = STR(int(IT),N+1)
      TF2 = STR(int(IT),N+2)
      TF3 = STR(int(IT),N+3)
      TF4 = STR(int(IT),N+4)
      TF5 = STR(int(IT),N+5)

      FUNCT = TF0 - DELT * (TF1 - DELT2 * (TF2 - DELT3 * (TF3 - DELT4 * &
                           (TF4 - DELT5 * TF5))))
   else
      FUNCT = FAC(N) * 1.D0 / (T**N * dsqrt(T))
   endif

end function funct
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
