!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! transform generic procedure
!
! HEADER DESCRIPTION PENDING
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
function transform_gen( Bmat, Cmat ) result( Dmat )
   implicit none
   GEN_TYPE        , intent(in)  :: Bmat(:,:)
   double precision, intent(in)  :: Cmat(:,:)
   GEN_TYPE        , allocatable :: Dmat(:,:), Xmat(:,:)
   double precision, allocatable :: Amat(:,:)
   logical :: error_found

   if (allocated(Amat)) deallocate(Amat)
   allocate(Amat(size(Cmat,2), size(Cmat,1)))
   Amat = transpose(Cmat)

#  include "matmul3_body.f90"
   ! This does not seem to work as intended. I am leaving it here for future
   ! reference -- F. Pedron Oct/2018

   !if (allocated(Mato)) deallocate(Mato)
   !allocate(Mato( size(Mati,1), size(Mati,2) ))
   !do jj = 1, size(Mati,2)
   !do ii = 1, size(Mati,1)
   !   Mato(ii,jj) = 0.0d0
   !   do kj = 1, size(Mati,2)
   !   do ki = 1, size(Mati,1)
   !      Mato(ii,jj) = Mato(ii,jj) + Bmat(ki,ii) * Mati(ki,kj) * Bmat(kj,jj)
   !   end do
   !   end do
   !end do
   !end do
   deallocate(Amat)
end function transform_gen
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
