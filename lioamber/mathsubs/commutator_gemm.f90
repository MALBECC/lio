!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!  CONMUTATOR - gemm version -
! 
!  CONMUTATOR(MA,MB)=MC=[MA*MB-MB*MA]
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function commutator_dd_gemm(MA,MB)
     > result(MC)
       implicit none
       real*8,intent(in)      :: MA(:,:)
       real*8,intent(in)      :: MB(:,:)
       real*8,allocatable     :: MC(:,:)
       real*8                 :: alpha,beta
       integer                :: nn
!
       nn=size(MA,1)
       allocate(MC(nn,nn))
       alpha=1.0D0
       beta=0.0D0
       call DGEMM('N','N',nn,nn,nn,alpha,MB,nn,MA,nn,beta,MC,nn)
       beta=-1.0D0
       call DGEMM('N','N',nn,nn,nn,alpha,MA,nn,MB,nn,beta,MC,nn)
       return;end function
!
!--------------------------------------------------------------------!
       function commutator_zd_gemm(MA,MB)
     > result(MC)
       implicit none
       complex*16,intent(in)  :: MA(:,:)
       real*8,intent(in)      :: MB(:,:)
       complex*16,allocatable :: MC(:,:)
       complex*16,allocatable :: MP(:,:)
       complex*16 :: alpha,beta
       integer                :: nn,i,j

       nn=size(MA,1)
       allocate(MC(nn,nn),MP(nn,nn))
       DO i=1,nn
          DO j=1,nn
             MP(i,j)=CMPLX(MB(i,j),0.0D0)
          ENDDO
       ENDDO
       alpha=CMPLX(1.0D0,0.0D0)
       beta=CMPLX(0.0D0,0.0D0)
!
       call ZGEMM('N','N',nn,nn,nn,alpha,MP,nn,MA,nn,beta,MC,nn)
       beta=CMPLX(-1.0D0,0.0D0)
       call ZGEMM('N','N',nn,nn,nn,alpha,MA,nn,MP,nn,beta,MC,nn)
!
       deallocate(MP)
       return;end function
!
!--------------------------------------------------------------------!
       function commutator_dz_gemm(MA,MB)
     > result(MC)
       implicit none
       real*8,intent(in)      :: MA(:,:)
       complex*16,intent(in)  :: MB(:,:)
       complex*16,allocatable :: MC(:,:)
       complex*16,allocatable :: MP(:,:)
       complex*16 :: alpha,beta
       integer                :: nn,i,j
!       
       nn=size(MA,1)
       allocate(MC(nn,nn),MP(nn,nn))
       DO i=1,nn
          DO j=1,nn
             MP(i,j)=CMPLX(MA(i,j),0.0D0)
          ENDDO
       ENDDO
       alpha=CMPLX(1.0D0,0.0D0)
       beta=CMPLX(0.0D0,0.0D0)
!      
       call ZGEMM('N','N',nn,nn,nn,alpha,MB,nn,MP,nn,beta,MC,nn)
       beta=CMPLX(-1.0D0,0.0D0)
       call ZGEMM('N','N',nn,nn,nn,alpha,MP,nn,MB,nn,beta,MC,nn) 
!
       deallocate(MP)      
       return;end function
!
!--------------------------------------------------------------------!
       function commutator_zz_gemm(MA,MB)
     > result(MC)
       implicit none
       complex*16,intent(in)  :: MA(:,:)
       complex*16,intent(in)  :: MB(:,:)
       complex*16,allocatable :: MC(:,:)
       complex*16 :: alpha,beta
       integer                :: nn

       nn=size(MA,1)
       allocate(MC(nn,nn))
       alpha=CMPLX(1.0D0,0.0D0)
       beta=CMPLX(0.0D0,0.0D0)
!      
       call ZGEMM('N','N',nn,nn,nn,alpha,MB,nn,MA,nn,beta,MC,nn)
       beta=CMPLX(-1.0D0,0.0D0)
       call ZGEMM('N','N',nn,nn,nn,alpha,MA,nn,MB,nn,beta,MC,nn)
!
       return;end function
!
!--------------------------------------------------------------------!
       function commutator_cd_gemm(MA,MB)
     > result(MC)
       implicit none
       complex*8,intent(in)  :: MA(:,:)
       real*8,intent(in)      :: MB(:,:)
       complex*8,allocatable :: MC(:,:)
       complex*8,allocatable :: MP(:,:)
       complex*8 :: alpha,beta
       integer                :: nn,i,j

       nn=size(MA,1)
       allocate(MC(nn,nn),MP(nn,nn))
       DO i=1,nn
          DO j=1,nn
             MP(i,j)=CMPLX(MA(i,j),0.0D0)
          ENDDO
       ENDDO
       alpha=CMPLX(1.0D0,0.0D0)
       beta=CMPLX(0.0D0,0.0D0)
!      
       call CGEMM('N','N',nn,nn,nn,alpha,MB,nn,MP,nn,beta,MC,nn)
       beta=CMPLX(-1.0D0,0.0D0)
       call CGEMM('N','N',nn,nn,nn,alpha,MP,nn,MB,nn,beta,MC,nn)
!
       deallocate(MP)
       return;end function
!
!--------------------------------------------------------------------!
       function commutator_dc_gemm(MA,MB)
     > result(MC)
       implicit none
       real*8,intent(in)      :: MA(:,:)
       complex*8,intent(in)  :: MB(:,:)
       complex*8,allocatable :: MC(:,:)
       complex*8,allocatable :: MP(:,:)
       complex*8 :: alpha,beta
       integer                :: nn,i,j
!
       nn=size(MA,1)
       allocate(MC(nn,nn),MP(nn,nn))
       DO i=1,nn
          DO j=1,nn
             MP(i,j)=CMPLX(MA(i,j),0.0D0)
          ENDDO
       ENDDO
       alpha=CMPLX(1.0D0,0.0D0)
       beta=CMPLX(0.0D0,0.0D0)    
       call CGEMM('N','N',nn,nn,nn,alpha,MB,nn,MP,nn,beta,MC,nn)
       beta=CMPLX(-1.0D0,0.0D0)
       call CGEMM('N','N',nn,nn,nn,alpha,MP,nn,MB,nn,beta,MC,nn)
       deallocate(MP)
       return;end function
!
!--------------------------------------------------------------------!
       function commutator_cc_gemm(MA,MB)
     > result(MC)
       implicit none
       complex*8,intent(in)  :: MA(:,:)
       complex*8,intent(in)  :: MB(:,:)
       complex*8,allocatable :: MC(:,:)
       complex*8 :: alpha,beta
       integer                :: nn

       nn=size(MA,1)
       allocate(MC(nn,nn))
       alpha=CMPLX(1.0D0,0.0D0)
       beta=CMPLX(0.0D0,0.0D0)
!      
       call CGEMM('N','N',nn,nn,nn,alpha,MB,nn,MA,nn,beta,MC,nn)
       beta=CMPLX(-1.0D0,0.0D0)
       call CGEMM('N','N',nn,nn,nn,alpha,MA,nn,MB,nn,beta,MC,nn)
       return;end function
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
