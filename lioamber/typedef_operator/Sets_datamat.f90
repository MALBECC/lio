subroutine Sets_data_AO (this, Dmat)

   implicit none
   class(operator), intent(inout) :: this
   real*8, intent(in)             :: Dmat(:,:)
   integer    :: Nbasis

   Nbasis = size( Dmat, 1 )

   if (allocated(this%data_AO)) then
      if (size(this%data_AO,1)/=Nbasis) then
        deallocate(this%data_AO)
        allocate(this%data_AO(Nbasis,Nbasis))
      endif
    else
      allocate(this%data_AO(Nbasis,Nbasis))
   endif

   this%data_AO = Dmat

end subroutine Sets_data_AO

subroutine Sets_data_ON (this, Dmat)

   implicit none
   class(operator), intent(inout) :: this
   real*8, intent(in)             :: Dmat(:,:)
   integer    :: Nbasis

   Nbasis = size( Dmat, 1 )

   if (allocated(this%data_ON)) then
      if (size(this%data_ON,1)/=Nbasis) then
        deallocate(this%data_ON)
        allocate(this%data_ON(Nbasis,Nbasis))
      endif
    else
      allocate(this%data_ON(Nbasis,Nbasis))
   endif

   this%data_ON = Dmat

end subroutine Sets_data_ON
