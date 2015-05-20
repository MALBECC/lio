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
       function basechange_d_cublas(M,Mati,devPtrX) result(Mato)
       implicit none
       integer, intent(in) :: M
       integer , intent(in)  :: devPtrX,devPtrXtrans
       REAL*8, intent(in) :: Mati(M,M)
       REAL*8, intent(out) :: Mato(M,M)
       REAL*8, dimension (:,:), ALLOCATABLE :: scratch
       allocate (scratch(M,M))
       call cumxtf(Mati,devPtrX,scratch,M)
       call cumfx(scratch,devPtrX,Mato,M)
       deallocate(scratch)
!
       return;end function
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function basechange_c_cublas(M,Mati,devPtrX,bctype) result(Mato)
       implicit none
       CHARACTER(LEN=3)   :: bctype
       integer*8 , intent(in)  :: devPtrX
       integer, intent(in) :: M
       COMPLEX*8, intent(in) :: Mati(M,M)
       COMPLEX*8, intent(out) :: Mato(M,M)
       COMPLEX*8, dimension (:,:), ALLOCATABLE :: scratch
       allocate (scratch(M,M))
       if(bctype.eq.'inv') then
          call cumxp(Mati,devPtrX,scratch,M)
          call cumpxt(scratch,devPtrX,Mato,M)
       elseif(bctype.eq.'dir') then
          call cumxtp(Mati,devPtrX,scratch,M)
          call cumpx(scratch,devPtrX,Mato,M)
       endif
       deallocate (scratch)
!
       return;end function
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function basechange_z_cublas(M,Mati,devPtrX,bctype) result(Mato)
       implicit none
       integer*8 , intent(in)  :: devPtrX
       CHARACTER(LEN=3)   :: bctype
       integer, intent(in) :: M
       COMPLEX*16, intent(in) :: Mati(M,M)
       COMPLEX*16, intent(out) :: Mato(M,M)
       COMPLEX*16, dimension (:,:), ALLOCATABLE :: scratch
       allocate (scratch(M,M))
       if(bctype.eq.'inv') then
          call cumxp(Mati,devPtrX,scratch,M)
          call cumpxt(scratch,devPtrX,Mato,M)
       elseif(bctype.eq.'dir') then
          call cumxtp(Mati,devPtrX,scratch,M)
          call cumpx(scratch,devPtrX,Mato,M)
       endif       
       deallocate (scratch)
!
       return;end function
!====================================================================!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
