!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% FUNC.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Dario Original - 11/March/1992                                               !
! This small program is taken from old version of gdfmol. It is the            !
! implementation of the Obara-Saika method for the evaluation of F(m,T), using !
! a 2 branch calculation. Ref: JCP 84 3963 (1986)                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

real*8 function funct(N, T)                                               
   use garcha_mod, only : STR, FAC

   implicit none
   real*8 , intent(inout) :: T
   integer, intent(in)    :: N

   integer :: it
   real*8  :: delt, delt2, delt3, delt4, delt5, fcap, ti, tf0, tf1, tf2, tf3, &
              tf4, tf5
     
   if (T.lt.0.0D0) then
      write(*,*) 'funct: Inout variable T is negative. Taking absolute value...'
      T=abs(T)
   endif
                                                                       
   if (T.le.43.975D0) then                          
      it = 20.0D0 * (T + 0.025D0)                                        
      ti = DFLOAT(it)                                                   
      it = it + 1                                                         

      delt  = T - 0.05D0 * ti                                             
      delt3 = delt   * 0.333333333333333D0                                
      delt4 = 0.25D0 * delt                                             
      delt2 = delt4  + delt4                                             
      delt5 = 0.20D0 * delt                                             

      tf0 = STR(it, N   )                                                    
      tf1 = STR(it, N +1)                                                    
      tf2 = STR(it, N +2)                                                    
      tf3 = STR(it, N +3)                                                    
      tf4 = STR(it, N +4)                                                    
      tf5 = STR(it, N +5)                                                    

      fcap  = tf0 -delt*(tf1 -delt2*(tf2 -delt3*(tf3 -delt4*(tf4 -delt5*tf5))))
      funct = fcap                                                   
   else   
      funct = FAC(N)*1.D0 / (T**N*dsqrt(T))
   endif  
   return
end function funct            
  
! Subroutine for table generation, used later in the Taylor expansion for 
! incomplete Gamma functions
subroutine generf
   use garcha_mod, only : STR, FAC
 
   implicit none
   real*8, parameter :: SQPI=1.77245385090551588D0

   real*8 :: t, u, w, y, fmch, f
   integer :: i, n, m

   ! Loops over t values in the table ( 0. to 43.95 , interval 0.05)
   do i = 1, 880
      t = 0.05D0 * DFLOAT(i-1)
      y = DEXP(-t)
      u = t + t
      f = fmch(22, t)
      ! Loop over order of incomple Gamma functions ( 0 to 21, the ones
      ! necessary for evaluating orders 0-16)
      do m = 21, 0, -1
         w = 2.D0* DFLOAT(m) + 1.D0
         str(i, m) = (y + u*f) /w
         f = STR(i, m)
      enddo
  enddo

  ! Calculation of the function [(v+1/2)/2**v+1*sqrt(pi) ]
  FAC(0) = SQPI /2.D0
  do n = 1, 16
      FAC(n) = FAC(n-1) * (2*n-1) /2
  enddo

  return
end subroutine generf

real*8 function fmch(M, X)
   implicit none
   real*8 , intent(in) :: X
   integer, intent(in) :: M

   real*8, parameter :: SQPI=1.77245385090551588D0
     
   real*8  ::  y, a, b, approx, fimult, fiprop, ptlsum, term, xd
   integer :: i, iw, notrms

   y = DEXP(-X)
   if (X.GT.10.D0) then
      a      = DFLOAT(M)
      b      = a + 0.5D0 
      xd     = 1.0D0 / X
      approx = 0.886226925452758D0 * DSQRT(xd) * xd**M

      do i = 1 , M
         approx = approx * (b - DFLOAT(i))
      enddo
  
      fimult = 0.5D0 * Y * xd
      ptlsum = 0.0D0

      if (fimult.eq.0.0D0) then 
         fmch = approx - fimult * ptlsum
         return
      endif

      fiprop = fimult / approx
      term   = 1.0D0
      ptlsum = term
      notrms = X + M

      do i = 2, notrms
         a      = b - DFLOAT(i-1)
         term   = term * a * xd
         ptlsum = ptlsum + term
         if( DABS(term*fiprop/ptlsum).le.1.0D-10 ) then
            fmch = approx - fimult * ptlsum
            return
         endif
      enddo
      
      write (iw, 999) M, X
      fmch = approx - fimult * ptlsum
      return

   else
      a      = DFLOAT(M) + 0.5D0
      term   = 1.0D0 / a
      ptlsum = term

      do i = 2, 50
         term = term * X / (A + DFLOAT(I-1))
         ptlsum = ptlsum + term
         if ( (term/ptlsum).lt.1.0D-12) then
            fmch = 0.5D0 * ptlsum * Y
            return
         endif
      enddo     
   endif

999 format ('fmch: convergence not achieved.', I6, 1PE16.8)
end function fmch
