!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function commutate_dd(MA,MB)
     > result(MC)
       implicit none
       real*8,intent(in)      :: MA(:,:)
       real*8,intent(in)      :: MB(:,:)
       real*8,allocatable     :: MC(:,:)
       real*8,allocatable     :: MP(:,:),MN(:,:)
       integer                :: nn

       nn=size(MA,1)
       allocate(MC(nn,nn),MP(nn,nn),MN(nn,nn))
       MP=MATMUL(MA,MB)
       MN=MATMUL(MB,MA)
       MC=MP-MN
       return;end function
!
!
!--------------------------------------------------------------------!
       function commutate_zd(MA,MB)
     > result(MC)
       implicit none
       complex*16,intent(in)  :: MA(:,:)
       real*8,intent(in)      :: MB(:,:)
       complex*16,allocatable :: MC(:,:)
       complex*16,allocatable :: MP(:,:),MN(:,:)
       integer                :: nn

       nn=size(MA,1)
       allocate(MC(nn,nn),MP(nn,nn),MN(nn,nn))
       MP=MATMUL(MA,MB)
       MN=MATMUL(MB,MA)
       MC=MP-MN
       return;end function
!
!
!--------------------------------------------------------------------!
       function commutate_dz(MA,MB)
     > result(MC)
       implicit none
       real*8,intent(in)      :: MA(:,:)
       complex*16,intent(in)  :: MB(:,:)
       complex*16,allocatable :: MC(:,:)
       complex*16,allocatable :: MP(:,:),MN(:,:)
       integer                :: nn

       nn=size(MA,1)
       allocate(MC(nn,nn),MP(nn,nn),MN(nn,nn))
       MP=MATMUL(MA,MB)
       MN=MATMUL(MB,MA)
       MC=MP-MN
       return;end function
!
!
!--------------------------------------------------------------------!
       function commutate_zz(MA,MB)
     > result(MC)
       implicit none
       complex*16,intent(in)  :: MA(:,:)
       complex*16,intent(in)  :: MB(:,:)
       complex*16,allocatable :: MC(:,:)
       complex*16,allocatable :: MP(:,:),MN(:,:)
       integer                :: nn

       nn=size(MA,1)
       allocate(MC(nn,nn),MP(nn,nn),MN(nn,nn))
       MP=MATMUL(MA,MB)
       MN=MATMUL(MB,MA)
       MC=MP-MN
       return;end function
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
