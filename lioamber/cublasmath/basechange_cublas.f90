!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! BASECHANGE - cublas version -
! BASETRANSFORM PROCEDURES
!
! (1) Initialization of scratch
! (2) First Product Xtrans*Mati
! (3) Second Product scratch*X
!
!====================================================================!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
function basechange_cublas_d(M, Mati, devPtrX, bctype) result(Mato)
   implicit none
   integer        , intent(in) :: M
   CUDAPTR, intent(in) :: devPtrX
   LIODBLE   , intent(in) :: Mati(M, M)

   character(len=3) :: bctype
   LIODBLE   , allocatable :: Mato(:,:), scratch(:,:)

   allocate (scratch(M, M), Mato(M, M))
   if (bctype == 'dir') then
      call cumxtf(Mati, devPtrX, scratch, M)
      call cumfx(scratch, devPtrX, Mato, M)
   elseif (bctype == 'inv') then
      call cumxtf(Mati, devPtrX, scratch, M)
      call cumfx(scratch, devPtrX, Mato, M)
   endif
   deallocate(scratch)
   
   return
end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
function basechange_cublas_c(M, Mati, devPtrX, bctype) result(Mato)
   implicit none
   integer        , intent(in) :: M 
   CUDAPTR, intent(in) :: devPtrX
   complex(kind=4), intent(in) :: Mati(M, M)
  
   character(len=3) :: bctype
   complex(kind=4), allocatable :: Mato(:,:), scratch(:,:)

   allocate (scratch(M, M), Mato(M, M))
   if (bctype == 'inv') then
      call cumxp(Mati, devPtrX, scratch, M)
      call cumpxt(scratch, devPtrX, Mato, M)
   elseif (bctype == 'dir') then
      call cumxtp(Mati, devPtrX, scratch, M)
      call cumpx(scratch, devPtrX, Mato, M)
   endif
   deallocate (scratch)

   return
end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
function basechange_cublas_z(M, Mati, devPtrX, bctype) result(Mato)
   implicit none
   integer        , intent(in) :: M
   CUDAPTR, intent(in) :: devPtrX
   complex(kind=8), intent(in) :: Mati(M, M)

   character(len=3) :: bctype
   complex(kind=8), allocatable :: Mato(:,:), scratch(:,:)

   allocate (scratch(M, M), Mato(M, M))
   if (bctype == 'inv') then
      call cumxp(Mati, devPtrX, scratch, M)
      call cumpxt(scratch, devPtrX, Mato, M)
   elseif (bctype == 'dir') then
      call cumxtp(Mati, devPtrX, scratch, M)
      call cumpx(scratch, devPtrX, Mato, M)
   endif   
   deallocate (scratch)
 
   return
end function
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
