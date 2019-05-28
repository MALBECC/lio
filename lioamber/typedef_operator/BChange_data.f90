!carlos: these subroutines allow the base changes of the matrix stored in
!        operator
subroutine BChange_AOtoON_r(this, Xmat, Nsize)
   use typedef_cumat, only: cumat_r

   implicit none
   class(operator) , intent(inout) :: this
   integer         , intent(in)    :: Nsize
   type(cumat_r)   , intent(in)    :: Xmat

   real(kind=8), allocatable :: Dmat(:,:)
   allocate(Dmat(Nsize,Nsize))
   Dmat = this%data_AO

   call Xmat%change_base(Dmat, 'dir')

   call this%Sets_data_ON(Dmat)
   deallocate(Dmat)
end subroutine BChange_AOtoON_r

subroutine BChange_ONtoAO_r(this, Xmat, Nsize)
   use typedef_cumat, only: cumat_r

   implicit none
   class(operator) , intent(inout) :: this
   integer         , intent(in)    :: Nsize
   type(cumat_r)   , intent(in)    :: Xmat

   real(kind=8), allocatable :: Dmat(:,:)
   
      allocate(Dmat(Nsize,Nsize))
      Dmat = this%data_ON

      call Xmat%change_base(Dmat, 'inv')

      call this%Sets_data_AO(Dmat)
      deallocate(Dmat)
 end subroutine BChange_ONtoAO_r

 subroutine BChange_AOtoON_x(this, Xmat, Nsize)
   use typedef_cumat, only: cumat_x

   implicit none
   class(operator) , intent(inout) :: this
   integer         , intent(in)    :: Nsize
   type(cumat_x)   , intent(in)    :: Xmat


   TDCOMPLEX   , allocatable :: DmatC(:,:)
   allocate(DmatC(Nsize,Nsize))
   DmatC = this%dataC_AO

   call Xmat%change_base(DmatC, 'dir')

   call this%Sets_dataC_ON(DmatC)
   deallocate(DmatC)

end subroutine BChange_AOtoON_x

subroutine BChange_ONtoAO_x(this, Xmat, Nsize)
   use typedef_cumat, only: cumat_x

   implicit none
   class(operator) , intent(inout) :: this
   integer         , intent(in)    :: Nsize
   type(cumat_x)   , intent(in)    :: Xmat

   TDCOMPLEX   , allocatable :: DmatC(:,:)
   
   
   allocate(DmatC(Nsize,Nsize))
   DmatC = this%dataC_ON
   call Xmat%change_base(DmatC, 'inv')

   call this%Sets_dataC_AO(DmatC)
   deallocate(DmatC)
   
 end subroutine BChange_ONtoAO_x