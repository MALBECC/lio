!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Sets_smat( this, Smat )
   
   implicit none
   class(sop), intent(inout) :: this
   real*8    , intent(in)    :: Smat(:,:)

   this%Nbasis = size( Smat, 1 )

   if ( allocated(this%Smat) ) deallocate(this%Smat)
   allocate( this%Smat( this%Nbasis, this%Nbasis ) )

   this%Smat = Smat

   if ( allocated(this%Lmat) ) deallocate(this%Lmat)
   allocate( this%Lmat( this%Nbasis, this%Nbasis ) )

   if ( allocated(this%Umat) ) deallocate(this%Umat)
   allocate( this%Umat( this%Nbasis, this%Nbasis ) )
   if ( allocated(this%Utrp) ) deallocate(this%Utrp)
   allocate( this%Utrp( this%Nbasis, this%Nbasis ) )

   if ( allocated(this%Gmat) ) deallocate(this%Gmat)
   allocate( this%Gmat( this%Nbasis, this%Nbasis ) )
   if ( allocated(this%Ginv) ) deallocate(this%Ginv)
   allocate( this%Ginv( this%Nbasis, this%Nbasis ) )

   if ( allocated(this%Vmat) ) deallocate(this%Vmat)
   allocate( this%Vmat( this%Nbasis, this%Nbasis ) )
   if ( allocated(this%Vtrp) ) deallocate(this%Vtrp)
   allocate( this%Vtrp( this%Nbasis, this%Nbasis ) )

   call this%Calc_fulldcmp_7m( this%Smat, this%Lmat, this%Umat, this%Utrp,     &
                             & this%Gmat, this%Ginv, this%Vmat, this%Vtrp )

end subroutine Sets_smat
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
