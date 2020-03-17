!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
function matmul3_ddd( Amat, Bmat, Cmat ) result( Dmat )
   implicit none
   LIODBLE, intent(in)  :: Amat(:,:)
   LIODBLE, intent(in)  :: Bmat(:,:)
   LIODBLE, intent(in)  :: Cmat(:,:)
   LIODBLE, allocatable :: Dmat(:,:)
   LIODBLE, allocatable :: Xmat(:,:)
   logical :: error_found
#  include "matmul3_body.f90"
end function matmul3_ddd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
function matmul3_dcd( A_in, Bmat, C_in ) result( Dmat )
   implicit none
   LIODBLE, intent(in)  :: A_in(:,:)
   TDCOMPLEX   , intent(in)  :: Bmat(:,:)
   LIODBLE, intent(in)  :: C_in(:,:)
   TDCOMPLEX, allocatable :: Dmat(:,:)
   TDCOMPLEX, allocatable :: Xmat(:,:), Amat(:,:), Cmat(:,:)
   logical :: error_found
   integer :: ii, jj
   TDCOMPLEX :: liocmplx

   ! This is necessary to avoid wrong type conversions in matmul.
   allocate(Amat(size(A_in,1), size(A_in,2)))
   do ii = 1, size(A_in,1)
   do jj = 1, size(A_in,2)
      Amat(ii,jj) = liocmplx(A_in(ii,jj),0.0D0)
   enddo
   enddo

   allocate(Amat(size(C_in,1), size(C_in,2)))
   do ii = 1, size(C_in,1)
   do jj = 1, size(C_in,2)
      Cmat(ii,jj) = liocmplx(C_in(ii,jj),0.0D0)
   enddo
   enddo

#  include "matmul3_body.f90"
   deallocate(Amat,Cmat)
end function matmul3_dcd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
