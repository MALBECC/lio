subroutine GenerateDensities(X,C,P,T,M,Ndim,NCO,Nvirt)
   implicit none

   integer, intent(in) :: M, Ndim, NCO, Nvirt
   double precision, intent(in)  :: X(Ndim), C(M,M)
   double precision, intent(out) :: P(M,M), T(M,M)

   double precision, dimension(:,:), allocatable :: PMO

   ! CALCULATE DIFFERENCE UNRELAXED DENSITY MATRIX
   allocate(PMO(M,M))
   call UnDiffDens(X,PMO,NCO,Nvirt,M,Ndim)
   call matMOtomatAO(PMO,P,C,M,M,.false.) ! TODO
   deallocate(PMO)

   ! CALCULATE TRANSITION DENSITY MATRIX
   call XmatForm(X,C,T,Ndim,NCO,Nvirt,M)
end subroutine GenerateDensities

subroutine UnDiffDens(X,T,NCO,Nvirt,M,Ndim)
   implicit none

   integer, intent(in) :: NCO, Nvirt, M, Ndim
   double precision, intent(in) :: X(Ndim)
   double precision, intent(out) :: T(M,M)

   integer :: i, j, pos, NCOc
   double precision, dimension(:,:), allocatable :: XM, XMtrans, Ptrash

   allocate(XM(NCO,Nvirt),XMtrans(Nvirt,NCO))

   T = 0.0D0
   NCOc = NCO + 1

   do i=1,NCO
   do j=1,Nvirt
     pos = (i-1) * Nvirt + j
     XM(i,j) = X(pos)
     XMtrans(j,i) = X(pos)
   enddo
   enddo

!  FORM UNRELAXED DIFFERENCE DENSITY MATRIX
   allocate(Ptrash(NCO,NCO))

!  FORM BLOCK OCC-OCC
   call dgemm('N','N',NCO,NCO,Nvirt,-1.0d0,XM,NCO,XMtrans, &
              Nvirt,0.0d0,Ptrash,NCO)

   do i=1,NCO
   do j=i,NCO
      T(NCOc-i,NCOc-j) = Ptrash(i,j)
      T(NCOc-j,NCOc-i) = Ptrash(i,j)
   enddo
   enddo

!  FORM BLOCK VIR - VIR
   deallocate(Ptrash); allocate(Ptrash(Nvirt,Nvirt))
   call dgemm('N','N',Nvirt,Nvirt,NCO,1.0d0,XMtrans,Nvirt,XM,&
              NCO,0.0d0,Ptrash,Nvirt)
 
   do i=1,Nvirt
   do j=i,Nvirt
      T(NCO+i,NCO+j) = Ptrash(i,j)
      T(NCO+j,NCO+i) = Ptrash(i,j)
   enddo
   enddo
   deallocate(Ptrash,XM,XMtrans)
end subroutine UnDiffDens

subroutine XmatForm(Vec,Coef,Mat,Ndim,NCO,Nvirt,M)
use excited_data, only: Coef_trans
   implicit none

   integer, intent(in) :: Ndim, NCO, Nvirt, M
   double precision, intent(in)  :: Vec(Ndim), Coef(M,M)
   double precision, intent(out) :: Mat(M,M)

   integer :: NCOc, row, col, pos
   double precision, dimension(:,:), allocatable ::SCR

   NCOc = NCO + 1
   do row=1,NCO
   do col=1,Nvirt
     pos = (row-1) * Nvirt + col
     Mat(NCOc-row,NCO+col) = Vec(pos)
   enddo
   enddo

   allocate(SCR(M,M))
   call dgemm('N','N',M,M,M,1.0d0,Coef,M,Mat,M,0.0d0,SCR,M)
   call dgemm('N','N',M,M,M,1.0d0,SCR,M,Coef_trans,M,0.0d0,Mat,M)
   deallocate(SCR)
end subroutine XmatForm
