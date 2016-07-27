!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! function matmult_XX( mata, matb, option, info ) result( matc )
!   implicit none
!   number(kind=X),intent(in)      :: mata(:,:)
!   number(kind=X),intent(in)      :: matb(:,:)
!   number(kind=X)                 :: matc(:,:)
!   integer, intent(in), optional  :: option
!   integer, intent(out), optional :: info

    matc = matmul( mata, matb )

! end function matmult_xx
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
