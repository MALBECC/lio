!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine calc_gintmat( mni, mnj, ai, aj, ri, rj, IntMat )
!------------------------------------------------------------------------------!
!
! IntMat(i,j) contains the integral of the product fa[i]*fb[j]
! between -infinity and +infinity, where:
!
!   fk[n](x) = (x-x[k])^n * exp( -alp[k] * (x-x[k])^2 )
!
! The matrix has an extra dimension of size 3 because there is
! one matrix for each of the 3 dimensions; the variables ai/aj
! are the same for all directions, but the x[k], y[k] and z[k]
! are different and are all introduced in ri[k],rj[k].
!
! TESTED
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in) :: mni,mnj
  real*8,intent(in)  :: ai,aj
  real*8,intent(in)  :: ri(3),rj(3)
  real*8,intent(out) :: IntMat(3,4,4)

  integer             :: ntot,ii,jj,kk
  real*8              :: a0,b1,b2,b3,b5,b7
  real*8              :: theta,PRE(3),r0(3)
  real*8              :: GI0,GI2,GI4,GI6
  real*8,dimension(3) :: di1,di2,di3,dj1,dj2,dj3
  real*8              :: term0,term2,term4

  real*8           :: SQPI
  real*8,parameter :: PI=3.14159265358979323846264338327d0



! Set initial parameters according to the matrix required
!------------------------------------------------------------------------------!
  SQPI=sqrt(PI)
  a0=ai+aj
  b2=(1/a0)

  b1=sqrt(b2)
  GI0=SQPI*b1

  b3=b1*b2
  GI2=SQPI*b3*(0.500d0)

  b5=b3*b2
  GI4=SQPI*b5*(0.750d0)

  b7=b5*b2
  GI6=SQPI*b7*(1.875d0)


  do kk=1,3
     r0(kk)=ri(kk)-rj(kk)
     r0(kk)=r0(kk)*ai*b2
     r0(kk)=r0(kk)+rj(kk)
     di1(kk)=r0(kk)-ri(kk)
     dj1(kk)=r0(kk)-rj(kk)

     theta=ri(kk)-rj(kk)
     theta=theta**2
     theta=theta*ai*aj
     theta=theta*b2
     PRE(kk)=exp(-theta)
  enddo

  do kk=1,3
     di2(kk)=di1(kk)*di1(kk)
     di3(kk)=di2(kk)*di1(kk)
     dj2(kk)=dj1(kk)*dj1(kk)
     dj3(kk)=dj2(kk)*dj1(kk)
  enddo


! Calculate the actual matrixes
!------------------------------------------------------------------------------!
  IntMat(:,:,:)=0.0d0
  do kk=1,3
     IntMat(kk,1,1)=GI0
     IntMat(kk,1,2)=GI0*dj1(kk)
     IntMat(kk,2,1)=GI0*di1(kk)
     IntMat(kk,2,2)=GI0*di1(kk)*dj1(kk)+GI2

     IntMat(kk,1,3)=GI0*dj2(kk)+GI2
     term0=GI0*di1(kk)*dj2(kk)
     term2=GI2*(di1(kk)+2*dj1(kk))
     IntMat(kk,2,3)=term0+term2

     IntMat(kk,3,1)=GI0*di2(kk)+GI2
     term0=GI0*di2(kk)*dj1(kk)
     term2=GI2*(2*di1(kk)+dj1(kk))
     IntMat(kk,3,2)=term0+term2

     term0=GI0*di2(kk)*dj2(kk)
     term2=GI2*(di2(kk)+4*di1(kk)*dj1(kk)+dj2(kk))
     IntMat(kk,3,3)=term0+term2+GI4

     term0=GI0*dj3(kk)
     term2=GI2*dj1(kk)*3
     IntMat(kk,1,4)=term0+term2

     term0=GI0*di1(kk)*dj3(kk)
     term2=GI2*3*(di1(kk)*dj1(kk)+dj2(kk))
     IntMat(kk,2,4)=term0+term2+GI4

     term0=GI0*di2(kk)*dj3(kk)
     term2=GI2*(3*di2(kk)*dj1(kk)+6*di1(kk)*dj2(kk)+dj3(kk))
     term4=GI4*(2*di1(kk)+3*dj1(kk))
     IntMat(kk,3,4)=term0+term2+term4

     term0=GI0*di3(kk)
     term2=GI2*di1(kk)*3
     IntMat(kk,4,1)=term0+term2

     term0=GI0*di3(kk)*dj1(kk)
     term2=GI2*3*(di2(kk)+di1(kk)*dj1(kk))
     IntMat(kk,4,2)=term0+term2+GI4

     term0=GI0*di3(kk)*dj2(kk)
     term2=GI2*(di3(kk)+6*di2(kk)*dj1(kk)+3*di1(kk)*dj2(kk))
     term4=GI4*(3*di1(kk)+2*dj1(kk))
     IntMat(kk,4,3)=term0+term2+term4

     term0=GI0*di3(kk)*dj3(kk)
     term2=di1(kk)*dj3(kk)+3*di2(kk)*dj2(kk)+di3(kk)*dj1(kk)
     term2=GI2*3*term2
     term4=GI4*3*(di2(kk)+3*di1(kk)*dj1(kk)+dj2(kk))
     IntMat(kk,4,4)=term0+term2+term4+GI6

     do jj=1,4
     do ii=1,4
        IntMat(kk,ii,jj)=IntMat(kk,ii,jj)*PRE(kk)
     enddo
     enddo


  enddo


end subroutine calc_gintmat
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
