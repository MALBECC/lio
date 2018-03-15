#ifdef CUBLAS
subroutine Commut_data(this, Bmat,devPtrX,devPtrY,AB_BAmat, Nsize)
   use cublasmath, only : cu_calc_fock_commuts
#else
subroutine Commut_data(this, Bmat,Xmat,Ymat,AB_BAmat, Nsize )
#endif

   implicit none
   class(operator), intent(inout) :: this
   integer, intent(in)            :: Nsize
   real*8, intent(in)             :: Bmat(Nsize,Nsize)

#ifdef  CUBLAS
!  ver intent pointers...
   integer*8, intent(in)         :: devPtrX
   integer*8, intent(in)         :: devPtrY
#else
   real*8,  intent(in)           :: Xmat(Nsize,Nsize)
   real*8,  intent(in)           :: Ymat(Nsize,Nsize)
#endif

   real*8, intent(out)           :: AB_BAmat(Nsize,Nsize)
   real*8   :: ABmat(Nsize,Nsize)
   real*8   :: BAmat(Nsize,Nsize)
   real*8   :: Amat(Nsize,Nsize)

   Amat=this%data_AO

#  ifdef  CUBLAS
      call cu_calc_fock_commuts(Amat,Bmat,devPtrX,devPtrY,AB_BAmat,Nsize)
#  else
      call calc_fock_commuts(Amat,Bmat,Xmat,Ymat,ABmat,BAmat,Nsize)
      AB_BAmat = ABmat - BAmat
#  endif

   call this%Sets_data_ON(Amat)

end subroutine Commut_data
