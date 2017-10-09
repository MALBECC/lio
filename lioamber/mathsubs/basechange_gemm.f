!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! BASECHANGE - manual version -
! BASETRANSFORM PROCEDURES
!
! (1) Initialization of Matm(nnd,ndd) and Mato(nii,ndd)
! (2) First Product Mati(nni,nnd)*Umat(nnd,ndd)
! (3) Second Product Utrp(nii,nni)*Matm(nni,ndd)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function basechange_d_gemm(M,Mati,Umat) result(Mato)
       implicit none
       integer,intent(in)     :: M
       real*8,intent(in)      :: Umat(M,M)
       real*8,intent(in)      :: Mati(M,M)
       real*8,allocatable     :: Matm(:,:)
       real*8,allocatable     :: Mato(:,:)
       integer                :: i,j
!
       allocate(Matm(M,M),Mato(M,M))
       Matm=DBLE(0)
       Mato=DBLE(0)
!
       call DGEMM('T','N',M,M,M,1.0D0,Umat,M,Mati,M,0.0D0,Matm,M)
       call DGEMM('N','N',M,M,M,1.0D0,Matm,M,Umat,M,0.0D0,Mato,M)
!
       return;end function
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function basechange_c_gemm(M,Mati,Umat) result(Mato)
       implicit none
       integer,intent(in)     :: M
       real*8,intent(in)      :: Umat(M,M)
       complex*8,allocatable  :: Cmat(:,:)
       complex*8,intent(in)   :: Mati(M,M)
       complex*8,allocatable  :: Matm(:,:)
       complex*8,allocatable  :: Mato(:,:)
       COMPLEX*8 :: alpha,beta
       integer                :: i,j
!
       allocate(Matm(M,M),Mato(M,M),Cmat(M,M))
       Matm=DCMPLX(0,0)
       Mato=DCMPLX(0,0)
       alpha=cmplx(1.0D0,0.0D0)
       beta=cmplx(0.0D0,0.0D0)
       DO i=1,M
          DO j=1,M
             Cmat(i,j)=cmplx(Umat(i,j),0.0D0)
          ENDDO
       ENDDO
!
       call CGEMM('T','N',M,M,M,alpha,Cmat,M,Mati,M,beta,Matm,M)
       call CGEMM('N','N',M,M,M,alpha,Matm,M,Cmat,M,beta,Mato,M)
!
       return;end function
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function basechange_z_gemm(M,Mati,Umat) result(Mato)
       implicit none
       integer,intent(in)      :: M
       real*8,intent(in)       :: Umat(M,M)
       complex*16,allocatable  :: Cmat(:,:)
       complex*16,intent(in)   :: Mati(M,M)
       complex*16,allocatable  :: Matm(:,:)
       complex*16,allocatable  :: Mato(:,:)
       complex*16              :: alpha,beta
       integer                 :: i,j
!
       allocate(Matm(M,M),Mato(M,M),Cmat(M,M))
       Matm=DCMPLX(0,0)
       Mato=DCMPLX(0,0)
       alpha=cmplx(1.0D0,0.0D0)
       beta=cmplx(0.0D0,0.0D0)
       DO i=1,M
          DO j=1,M
             Cmat(i,j)=cmplx(Umat(i,j),0.0D0)
          ENDDO
       ENDDO
!
       call ZGEMM('T','N',M,M,M,alpha,Cmat,M,Mati,M,beta,Matm,M)
       call ZGEMM('N','N',M,M,M,alpha,Matm,M,Cmat,M,beta,Mato,M)
!
       return;end function
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
