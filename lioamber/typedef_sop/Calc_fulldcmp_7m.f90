!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Calc_fulldcmp_7m( Smat, Lmat, Umat, Utrp, Gmat, Ginv, Vmat, Vtrp )
!------------------------------------------------------------------------------!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   implicit none
   LIODBLE , intent(in)  :: Smat(:,:)
   LIODBLE , intent(out) :: Lmat(:,:)
   LIODBLE , intent(out) :: Umat(:,:)
   LIODBLE , intent(out) :: Utrp(:,:)
   LIODBLE , intent(out) :: Gmat(:,:)
   LIODBLE , intent(out) :: Ginv(:,:)
   LIODBLE , intent(out) :: Vmat(:,:)
   LIODBLE , intent(out) :: Vtrp(:,:)

   integer              :: Nbasis, nn

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
