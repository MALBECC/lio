
#ifdef CUBLAS
   subroutine BChange_AOtoON(this,devPtrX,Nsize)
   use cublasmath, only : cumfx, cumxtf
#else
   subroutine BChange_AOtoON(this,Xmat,Nsize)
   use mathsubs, only : basechange_gemm
#endif

   implicit none
   class(operator), intent(inout) :: this
   integer, intent(in)           :: Nsize
#  ifdef  CUBLAS
   integer*8, intent(in) :: devPtrX
#  else
   real*8,  intent(in) :: Xmat(Nsize,Nsize)
#  endif

   real*8 :: Dmat(Nsize,Nsize)

   Dmat = this%data_AO

#  ifdef  CUBLAS
     call cumxtf(Dmat,devPtrX,Dmat,Nsize)
     call cumfx(Dmat,DevPtrX,Dmat,Nsize)
#  else
      Dmat = basechange_gemm(Nsize,Dmat,Xmat)
#  endif

   call this%Sets_data_ON(Dmat)

 end subroutine BChange_AOtoON
