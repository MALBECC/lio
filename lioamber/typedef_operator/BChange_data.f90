!carlos: these subroutines allow the base changes of the matrix stored in
!        operator
subroutine BChange_AOtoON_r(this, Xmat, Nsize, mode)
   use typedef_cumat, only: cumat_r

   implicit none
   class(operator) , intent(inout) :: this
   character(len=1), intent(in)    :: mode
   integer         , intent(in)    :: Nsize
   type(cumat_r)   , intent(in)    :: Xmat


   real(kind=8), allocatable :: Dmat(:,:)
   TDCOMPLEX   , allocatable :: DmatC(:,:)

   if (mode == 'r') then
      allocate(Dmat(Nsize,Nsize))
      Dmat = this%data_AO

      call Xmat%change_base(Dmat, 'dir')

      call this%Sets_data_ON(Dmat)
      deallocate(Dmat)

   elseif (mode == 'c') then
      allocate(DmatC(Nsize,Nsize))
      DmatC = this%dataC_AO

      call Xmat%change_base(DmatC, 'dir')

      call this%Sets_dataC_ON(DmatC)
      deallocate(DmatC)
   endif

end subroutine BChange_AOtoON_r

subroutine BChange_ONtoAO_r(this, Xmat, Nsize, mode)
   use typedef_cumat, only: cumat_r

   implicit none
   class(operator) , intent(inout) :: this
   character(len=1), intent(in)    :: mode
   integer         , intent(in)    :: Nsize
   type(cumat_r)   , intent(in)    :: Xmat

   real(kind=8), allocatable :: Dmat(:,:)
   TDCOMPLEX   , allocatable :: DmatC(:,:)
   
   if (mode == 'r') then
      allocate(Dmat(Nsize,Nsize))
      Dmat = this%data_AO

      call Xmat%change_base(Dmat, 'inv')

      call this%Sets_data_ON(Dmat)
      deallocate(Dmat)

   elseif (mode == 'c') then
      allocate(DmatC(Nsize,Nsize))
      DmatC = this%dataC_AO

      call Xmat%change_base(DmatC, 'inv')

      call this%Sets_dataC_ON(DmatC)
      deallocate(DmatC)
   endif

 end subroutine BChange_ONtoAO_r

 subroutine BChange_AOtoON_x(this, Xmat, Nsize, mode)
   use typedef_cumat, only: cumat_x

   implicit none
   class(operator) , intent(inout) :: this
   character(len=1), intent(in)    :: mode
   integer         , intent(in)    :: Nsize
   type(cumat_x)   , intent(in)    :: Xmat


   real(kind=8), allocatable :: Dmat(:,:)
   TDCOMPLEX   , allocatable :: DmatC(:,:)

   if (mode == 'r') then
      allocate(Dmat(Nsize,Nsize))
      Dmat = this%data_AO

      call Xmat%change_base(Dmat, 'dir')

      call this%Sets_data_ON(Dmat)
      deallocate(Dmat)

   elseif (mode == 'c') then
      allocate(DmatC(Nsize,Nsize))
      DmatC = this%dataC_AO

      call Xmat%change_base(DmatC, 'dir')

      call this%Sets_dataC_ON(DmatC)
      deallocate(DmatC)
   endif

end subroutine BChange_AOtoON_x

subroutine BChange_ONtoAO_x(this, Xmat, Nsize, mode)
   use typedef_cumat, only: cumat_x

   implicit none
   class(operator) , intent(inout) :: this
   character(len=1), intent(in)    :: mode
   integer         , intent(in)    :: Nsize
   type(cumat_x)   , intent(in)    :: Xmat

   real(kind=8), allocatable :: Dmat(:,:)
   TDCOMPLEX   , allocatable :: DmatC(:,:)
   
   if (mode == 'r') then
      allocate(Dmat(Nsize,Nsize))
      Dmat = this%data_AO

      call Xmat%change_base(Dmat, 'inv')

      call this%Sets_data_ON(Dmat)
      deallocate(Dmat)

   elseif (mode == 'c') then
      allocate(DmatC(Nsize,Nsize))
      DmatC = this%dataC_AO

      call Xmat%change_base(DmatC, 'inv')

      call this%Sets_dataC_ON(DmatC)
      deallocate(DmatC)
   endif

 end subroutine BChange_ONtoAO_x