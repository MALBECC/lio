!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! transform generic procedure
!
! HEADER DESCRIPTION PENDING
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
function transform_gen( Bmat, C_in ) result( Dmat )
   implicit none
   GEN_TYPE        , intent(in)  :: Bmat(:,:)
   LIODBLE, intent(in)  :: C_in(:,:)
   GEN_TYPE        , allocatable :: Dmat(:,:), Xmat(:,:)
   GEN_TYPE        , allocatable :: Amat(:,:), Cmat(:,:)
   logical :: error_found
#ifdef CONVERT_R
   integer :: ii, jj
#endif

   if (allocated(Amat)) deallocate(Amat)
   allocate(Amat(size(C_in,2), size(C_in,1)), Cmat(size(C_in,2), size(C_in,1)))

#ifdef CONVERT_R
   do ii = 1, size(C_in,1)
   do jj = 1, size(C_in,2)
       Cmat(ii,jj) = cmplx(C_in(ii,jj),0.0D0,CPSIZE)
   enddo
   enddo
#else
   Cmat = C_in
#endif

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
   deallocate(Amat, Cmat)
end function transform_gen
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
