!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function gaussbell_d(pos,mid,alpha,vector0)
     > result(vector)
       implicit none
       real*8,intent(in)     :: pos,mid,alpha
       real*8,intent(in)     :: vector0(:)
       real*8                :: vector(size(vector0))
       integer               :: kk
       real*8                :: factor

       factor=exp((-1)*(alpha)*((pos-mid)**2))
       do kk=1,size(vector0)
         vector(kk)=vector0(kk)*factor
       enddo
       return;end function
!
!
!--------------------------------------------------------------------!
       function gaussbell_z(pos,mid,alpha,vector0)
     > result(vector)
       implicit none
       real*8,intent(in)     :: pos,mid,alpha
       complex*16,intent(in) :: vector0(:)
       complex*16            :: vector(size(vector0))
       integer               :: kk
       real*8                :: factor

       factor=exp((-1)*(alpha)*((pos-mid)**2))
       do kk=1,size(vector0)
         vector(kk)=vector0(kk)*factor
       enddo
       return;end function
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
