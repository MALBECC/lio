!carlos: for the moment these subroutines only conmut the ON and real type matrix
subroutine Commut_data_r(this, Bmat, AB_BAmat, Nsize)

#ifdef CUBLAS
   use cublasmath, only : commutator_cublas
#else
   use mathsubs,   only: commutator
#endif

   implicit none
   class(operator), intent(inout) :: this
   integer, intent(in)            :: Nsize
   real*8, intent(in)             :: Bmat(Nsize,Nsize)
   real*8, intent(inout)           :: AB_BAmat(Nsize,Nsize)

   real*8, allocatable  :: ABmat(:,:)
   real*8, allocatable  :: BAmat(:,:)
   real*8, allocatable :: Amat(:,:)
   allocate(Amat(Nsize,Nsize), ABmat(Nsize,Nsize),BAmat(Nsize,Nsize))

   Amat=this%data_ON

#ifdef CUBLAS
      AB_BAmat = commutator_cublas(Amat, Bmat)
#else
      AB_BAmat = commutator (Amat, Bmat)
#endif

end subroutine Commut_data_r

subroutine Commut_data_c(this, Bmat, AB_BAmat, Nsize)

#ifdef CUBLAS
   use cublasmath, only : commutator_cublas
#else
   use mathsubs,   only: commutator
#endif

   implicit none
   class(operator), intent(inout) :: this
   integer, intent(in)            :: Nsize

   TDCOMPLEX, intent(in)       :: Bmat(Nsize,Nsize)
   TDCOMPLEX, intent(out)      :: AB_BAmat(Nsize,Nsize)
   TDCOMPLEX, allocatable      :: ABmat(:,:)
   TDCOMPLEX, allocatable      :: BAmat(:,:)

   real*8, allocatable :: Amat(:,:)

   allocate(Amat(Nsize,Nsize), ABmat(Nsize,Nsize), BAmat(Nsize,Nsize))

   Amat=this%data_ON

#ifdef CUBLAS
      AB_BAmat = commutator_cublas(Amat, Bmat)
#else
      AB_BAmat = commutator (Amat, Bmat)
#endif

end subroutine Commut_data_c
