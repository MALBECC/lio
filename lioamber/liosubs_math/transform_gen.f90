!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! transform generic procedure
!
! HEADER DESCRIPTION PENDING
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
function transform_gen( Mati, Bmat ) result( Mato )
   implicit none
   GEN_TYPE, intent(in) :: Mati(:,:)
   real*8, intent(in)   :: Bmat(:,:)
   GEN_TYPE             :: Mato( size(Mati,1), size(Mati,2) )
   integer              :: ii, jj, ki, kj

   do jj = 1, size(Mati,2)
   do ii = 1, size(Mati,1)
      Mato(ii,jj) = 0.0d0
      do kj = 1, size(Mati,2)
      do ki = 1, size(Mati,1)
         Mato(ii,jj) = Mato(ii,jj) + Bmat(ki,ii) * Mati(ki,kj) * Bmat(kj,jj)
      end do
      end do
   end do
   end do


end function transform_gen

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
