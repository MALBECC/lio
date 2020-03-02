!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Calc_fulldcmp_4m( Smat, Lmat, Umat, Gmat, Vtrp )
!------------------------------------------------------------------------------!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   use liosubs_math, only: matdcmp_cholesky, matdcmp_svd
   implicit none
   LIODBLE , intent(in)  :: Smat(:,:)
   LIODBLE , intent(out) :: Lmat(:,:)
   LIODBLE , intent(out) :: Umat(:,:)
   LIODBLE , intent(out) :: Gmat(:,:)
   LIODBLE , intent(out) :: Vtrp(:,:)

   integer              :: Nbasis

   Nbasis = size( Smat, 1 )
   ! TODO: size consistency checks

   call matdcmp_cholesky( Smat, Lmat )

   call matdcmp_svd( Lmat, Umat, Gmat, Vtrp )

end subroutine Calc_fulldcmp_4m
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
