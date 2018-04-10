!carlos: this subroutine builds the density matrix.
subroutine Dens_build(this,Msize, Nocup,Focup,coef_mat)
   use liosubs_dens  , only: builds_densmat

   implicit none
   class(operator), intent(inout) :: this
   integer, intent(in)            :: Msize
   integer, intent(in)            :: Nocup
   real*8 , intent(in)            :: Focup
   real*8 , intent(in)            :: coef_mat(Msize, Msize)

   call builds_densmat( Msize, Nocup, Focup, coef_mat, this%data_AO)

end subroutine Dens_build
