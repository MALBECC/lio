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
       function basechange_cublas_d(M,Mati,devPtrX,bctype) result(Mato)
       implicit none
       CHARACTER(LEN=3)   :: bctype
       integer, intent(in) :: M
       integer*8, intent(in)  :: devPtrX
       REAL*8, intent(in) :: Mati(M,M)
       REAL*8, ALLOCATABLE :: Mato(:,:)
       REAL*8, dimension (:,:), ALLOCATABLE :: scratch
       allocate (scratch(M,M),Mato(M,M))
       if(bctype.eq.'dir') then
       call cumxtf(Mati,devPtrX,scratch,M)
       call cumfx(scratch,devPtrX,Mato,M)
       elseif(bctype.eq.'inv') then
       call cumxtf(Mati,devPtrX,scratch,M)
       call cumfx(scratch,devPtrX,Mato,M)
       endif
       deallocate(scratch)
!
       return;end function
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function basechange_cublas_c(M,Mati,devPtrX,bctype) result(Mato)
       implicit none
       CHARACTER(LEN=3)   :: bctype
       integer*8 , intent(in)  :: devPtrX
       integer, intent(in) :: M
       COMPLEX*8, intent(in) :: Mati(M,M)
       COMPLEX*8, ALLOCATABLE :: Mato(:,:)
       COMPLEX*8, dimension (:,:), ALLOCATABLE :: scratch
       allocate (scratch(M,M),Mato(M,M))
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
       function basechange_cublas_z(M,Mati,devPtrX,bctype) result(Mato)
       implicit none
       integer*8 , intent(in)  :: devPtrX
       CHARACTER(LEN=3)   :: bctype
       integer, intent(in) :: M
       COMPLEX*16, intent(in) :: Mati(M,M)
       COMPLEX*16, ALLOCATABLE :: Mato(:,:)
       COMPLEX*16, dimension (:,:), ALLOCATABLE :: scratch
       allocate (scratch(M,M),Mato(M,M))
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
