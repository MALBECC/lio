!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! BASECHANGE - manual version -
! BASETRANSFORM PROCEDURES
!
! (1) Initialization of Matm(nnd,ndd) and Mato(nii,ndd)
! (2) First Product Mati(nni,nnd)*Umat(nnd,ndd)
! (3) Second Product Utrp(nii,nni)*Matm(nni,ndd)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function basechange_d_gemm(M,Mati,Umat,mode) result(Mato)
       implicit none
       integer,intent(in)     :: M
       character(len=3)       :: mode
       real*8,intent(in)      :: Umat(M,M)
       real*8,intent(in)      :: Mati(M,M)
       real*8,allocatable     :: Matm(:,:)
       real*8,allocatable     :: Mato(:,:)
       
       allocate(Matm(M,M), Mato(M,M))
       Matm = 0.0D0
       Mato = 0.0D0

       if (mode == 'inv') then
          call DGEMM('N','N',M,M,M,1.0D0,Umat,M,Mati,M,0.0D0,Matm,M)
          call DGEMM('N','T',M,M,M,1.0D0,Matm,M,Umat,M,0.0D0,Mato,M)
       else
          call DGEMM('T','N',M,M,M,1.0D0,Umat,M,Mati,M,0.0D0,Matm,M)
          call DGEMM('N','N',M,M,M,1.0D0,Matm,M,Umat,M,0.0D0,Mato,M)
       endif
       deallocate(Matm)

       return
       end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function basechange_cc_gemm(M,Mati,Umat,mode) result(Mato)
       implicit none
       integer,intent(in)     :: M
       character(len=3)       :: mode
       complex*8,intent(in)   :: Umat(M,M)
       complex*8,intent(in)   :: Mati(M,M)
       complex*8,allocatable  :: Matm(:,:)
       complex*8,allocatable  :: Mato(:,:)
       COMPLEX*8 :: alpha,beta

       allocate(Matm(M,M),Mato(M,M))
       Matm  = cmplx(0.0E0,0.0E0)
       Mato  = cmplx(0.0E0,0.0E0)
       alpha = cmplx(1.0E0,0.0E0)
       beta  = cmplx(0.0E0,0.0E0)
       if (mode == 'inv') then
          call CGEMM('N','N',M,M,M,alpha,Umat,M,Mati,M,beta,Matm,M)
          call CGEMM('N','T',M,M,M,alpha,Matm,M,Umat,M,beta,Mato,M)
       else
          call CGEMM('T','N',M,M,M,alpha,Umat,M,Mati,M,beta,Matm,M)
          call CGEMM('N','N',M,M,M,alpha,Matm,M,Umat,M,beta,Mato,M)
       endif
       deallocate(Matm)

       return
       end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function basechange_zz_gemm(M,Mati,Umat,mode) result(Mato)
       implicit none
       integer,intent(in)      :: M
       character(len=3)        :: mode
       complex*16,intent(in)   :: Umat(M,M)
       complex*16,intent(in)   :: Mati(M,M)
       complex*16,allocatable  :: Matm(:,:)
       complex*16,allocatable  :: Mato(:,:)
       complex*16              :: alpha,beta

       allocate(Matm(M,M),Mato(M,M))
       Matm  = dcmplx(0.0D0,0.0D0)
       Mato  = dcmplx(0.0D0,0.0D0)
       alpha = dcmplx(1.0D0,0.0D0)
       beta  = dcmplx(0.0D0,0.0D0)

       if (mode == 'inv') then
          call ZGEMM('N','N',M,M,M,alpha,Umat,M,Mati,M,beta,Matm,M)
          call ZGEMM('N','T',M,M,M,alpha,Matm,M,Umat,M,beta,Mato,M)
       else
          call ZGEMM('T','N',M,M,M,alpha,Umat,M,Mati,M,beta,Matm,M)
          call ZGEMM('N','N',M,M,M,alpha,Matm,M,Umat,M,beta,Mato,M)
       endif
       deallocate(Matm)

       return
       end function
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
