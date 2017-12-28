!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Calc_fulldcmp_7m( Smat, Lmat, Umat, Utrp, Gmat, Ginv, Vmat, Vtrp )
!------------------------------------------------------------------------------!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   implicit none
   real*8 , intent(in)  :: Smat(:,:)
   real*8 , intent(out) :: Lmat(:,:)
   real*8 , intent(out) :: Umat(:,:)
   real*8 , intent(out) :: Utrp(:,:)
   real*8 , intent(out) :: Gmat(:,:)
   real*8 , intent(out) :: Ginv(:,:)
   real*8 , intent(out) :: Vmat(:,:)
   real*8 , intent(out) :: Vtrp(:,:)

   integer              :: Nbasis, nn, ii, jj

   Nbasis = size( Smat, 1 )
   ! TODO: size consistency checks

   call Calc_fulldcmp_4m( Smat, Lmat, Umat, Gmat, Vtrp )

   Utrp = transpose( Umat )
   Vmat = transpose( Vtrp )

   Ginv(:,:) = 0.0d0
   do nn = 1, Nbasis
      Ginv(nn,nn) = 1 / Gmat(nn,nn)
   enddo

end subroutine Calc_fulldcmp_7m
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
