!carlos: these subroutines store matrices in operator.
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Sets_dataC_AO (this, Dmat)

   implicit none
   class(operator), intent(inout)  :: this
#ifdef TD_SIMPLE
   complex*8, intent(in)           :: Dmat(:,:)
#else
    complex*16, intent(in)         :: Dmat(:,:)
#endif
   integer    :: Nbasis

   Nbasis = size( Dmat, 1 )

   if (allocated(this%dataC_AO)) then
      if (size(this%dataC_AO,1)/=Nbasis) then
        deallocate(this%dataC_AO)
        allocate(this%dataC_AO(Nbasis,Nbasis))
      endif
    else
      allocate(this%dataC_AO(Nbasis,Nbasis))
   endif

   this%dataC_AO = Dmat

end subroutine Sets_dataC_AO
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Sets_dataC_ON (this, Dmat)

   implicit none
   class(operator), intent(inout)  :: this
#ifdef TD_SIMPLE
   complex*8, intent(in)           :: Dmat(:,:)
#else
    complex*16, intent(in)         :: Dmat(:,:)
#endif
   integer    :: Nbasis

   Nbasis = size( Dmat, 1 )

   if (allocated(this%dataC_ON)) then
      if (size(this%dataC_ON,1)/=Nbasis) then
        deallocate(this%dataC_ON)
        allocate(this%dataC_ON(Nbasis,Nbasis))
      endif
    else
      allocate(this%dataC_ON(Nbasis,Nbasis))
   endif

   this%dataC_ON = Dmat

end subroutine Sets_dataC_ON
