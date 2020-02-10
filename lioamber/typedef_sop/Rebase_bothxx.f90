!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Rebase_fockao( this, fock )
   use liosubs_math, only: matmul3
   implicit none
   class(sop), intent(in)    :: this
   real*8    , intent(inout) :: fock( this%Nbasis, this%Nbasis )

!  From AO to ON
   fock = matmul3( this%Xtrp, fock, this%Xmat )

end subroutine Rebase_fockao
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Rebase_fockon( this, fock )
   use liosubs_math, only: matmul3
   implicit none
   class(sop), intent(in)    :: this
   real*8    , intent(inout) :: fock( this%Nbasis, this%Nbasis )

!  From ON to AO
   fock = matmul3( this%Ymat, fock, this%Ytrp )

end subroutine Rebase_fockon
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Rebase_densao( this, dens )
   use liosubs_math, only: matmul3
   implicit none
   class(sop), intent(in)    :: this
   complex*16, intent(inout) :: dens( this%Nbasis, this%Nbasis )

!  From AO to ON
   dens = matmul3( this%Ytrp, dens, this%Ymat )

end subroutine Rebase_densao
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Rebase_denson( this, dens )
   use liosubs_math, only: matmul3
   implicit none
   class(sop), intent(in)    :: this
   complex*16, intent(inout) :: dens( this%Nbasis, this%Nbasis )

!  From ON to AO
   dens = matmul3( this%Xmat, dens, this%Xtrp )

end subroutine Rebase_denson
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
