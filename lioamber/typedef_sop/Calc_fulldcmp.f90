!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Calc_fulldcmp( Smat, Umat, Utrp, Gmat, Ginv, Vmat, Vtrp )
!------------------------------------------------------------------------------!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   use liosubs_math, only: matdcmp_cholesky, matdcmp_svd
   implicit none
   real*8 , intent(in)  :: Smat(:,:)
   real*8 , intent(out) :: Umat(:,:)
   real*8 , intent(out) :: Utrp(:,:)
   real*8 , intent(out) :: Gmat(:,:)
   real*8 , intent(out) :: Ginv(:,:)
   real*8 , intent(out) :: Vmat(:,:)
   real*8 , intent(out) :: Vtrp(:,:)

   integer              :: Nbasis, nn, ii, jj
   real*8 , allocatable :: Lmat(:,:)

   Nbasis = size( Smat, 1 )
   ! TODO: size consistency checks

   allocate( Lmat( Nbasis, Nbasis ) )

   call matdcmp_cholesky( Smat, Lmat )

   call matdcmp_svd( Lmat, Umat, Gmat, Vtrp )

   Ginv(:,:) = 0.0d0
   do nn = 1, Nbasis
      Ginv(nn,nn) = 1 / Gmat(nn,nn)
   enddo

   Utrp = transpose( Umat )
   Vmat = transpose( Vtrp )

   deallocate( Lmat )
end subroutine Calc_fulldcmp
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
