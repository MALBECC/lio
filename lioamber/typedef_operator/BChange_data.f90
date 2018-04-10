!carlos: these subroutines allow the base changes of the matrix stored in
!        operator
#ifdef CUBLAS
   subroutine BChange_AOtoON(this,devPtrX,Nsize,mode)
   use cublasmath, only : basechange_cublas
#else
   subroutine BChange_AOtoON(this,Xmat,Nsize, mode)
   use mathsubs, only : basechange_gemm
#endif

   implicit none
   class(operator), intent(inout) :: this
   character(len=1), intent(in)    :: mode
   integer, intent(in)            :: Nsize
#ifdef  CUBLAS
   integer*8, intent(in) :: devPtrX
#else
   real*8,  intent(in) :: Xmat(Nsize,Nsize)
#endif

   real*8, allocatable :: Dmat(:,:)
#ifdef TD_SIMPLE
   complex*8, allocatable  :: DmatC(:,:)
#else
   complex*16, allocatable :: DmatC(:,:)
#endif

   if (mode.eq.'r') then
      allocate(Dmat(Nsize,Nsize))
      Dmat = this%data_AO
#ifdef  CUBLAS
      Dmat=basechange_cublas(Nsize,Dmat,devPtrX,'dir')
#else
      Dmat = basechange_gemm(Nsize,Dmat,Xmat)
#endif

      call this%Sets_data_ON(Dmat)

   else if (mode.eq.'c') then
      allocate(DmatC(Nsize,Nsize))
      DmatC = this%dataC_AO
#ifdef  CUBLAS
      DmatC=basechange_cublas(Nsize,DmatC,devPtrX,'dir')
#else
      DmatC= basechange_gemm(Nsize,DmatC,Xmat)
#endif
      call this%Sets_dataC_ON(DmatC)
   end if

 end subroutine BChange_AOtoON

#ifdef CUBLAS
   subroutine BChange_ONtoAO(this,devPtrX,Nsize,mode)
   use cublasmath, only : basechange_cublas
#else
   subroutine BChange_ONtoAO(this,Xmat,Nsize, mode)
   use mathsubs, only : basechange_gemm
#endif

   implicit none
   class(operator), intent(inout) :: this
   character(len=1), intent(in)    :: mode
   integer, intent(in)            :: Nsize
#ifdef  CUBLAS
   integer*8, intent(in) :: devPtrX
#else
   real*8,  intent(in) :: Xmat(Nsize,Nsize)
#endif

   real*8, allocatable :: Dmat(:,:)
#ifdef TD_SIMPLE
   complex*8, allocatable  :: DmatC(:,:)
#else
   complex*16, allocatable :: DmatC(:,:)
#endif

   if (mode.eq.'r') then
      allocate(Dmat(Nsize,Nsize))
      Dmat = this%data_ON
#ifdef  CUBLAS
!charly: estas probando basechange cublas, despues fijate si anda que no compilaste
      Dmat=basechange_cublas(Nsize,Dmat,devPtrX,'inv')
#else
      Dmat = basechange_gemm(Nsize,Dmat,Xmat)
#endif

      call this%Sets_data_AO(Dmat)

   else if (mode.eq.'c') then
      allocate(DmatC(Nsize,Nsize))
      DmatC = this%dataC_ON
#ifdef  CUBLAS
      DmatC=basechange_cublas(Nsize,DmatC,devPtrX,'inv')
#else
      DmatC= basechange_gemm(Nsize,DmatC,Xmat)
#endif
      call this%Sets_dataC_AO(DmatC)
   end if

 end subroutine BChange_ONtoAO
