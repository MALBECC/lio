!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Calc_fulldcmp_4m( Smat, Lmat, Umat, Gmat, Vtrp )
!------------------------------------------------------------------------------!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   use liosubs_math, only: matdcmp_cholesky, matdcmp_svd
   implicit none
   real*8 , intent(in)  :: Smat(:,:)
   real*8 , intent(out) :: Lmat(:,:)
   real*8 , intent(out) :: Umat(:,:)
   real*8 , intent(out) :: Gmat(:,:)
   real*8 , intent(out) :: Vtrp(:,:)

   integer              :: Nbasis, nn, ii, jj

   Nbasis = size( Smat, 1 )
   ! TODO: size consistency checks

   call matdcmp_cholesky( Smat, Lmat )

   call matdcmp_svd( Lmat, Umat, Gmat, Vtrp )

end subroutine Calc_fulldcmp_4m
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
