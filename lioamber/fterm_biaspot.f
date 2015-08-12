!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine fterm_biaspot(M,sqsmat,vector,weight,outmat)
       implicit none
       integer,intent(in)   :: M
       real*8,intent(in)    :: sqsmat(M,M)
       integer,intent(in)   :: vector(M)
       real*8,intent(in)    :: weight
       real*8,intent(inout) :: outmat(M,M)

       real*8  :: newterm,foruri
       integer :: ii,jj,kk

       do ii=1,M
       do jj=1,M
!         foruri=0.0d0
         do kk=1,M
           newterm=sqsmat(ii,kk)*weight*sqsmat(kk,jj)
           outmat(ii,jj)=outmat(ii,jj)+newterm*real(vector(kk))
!           foruri=foruri+sqsmat(ii,kk)*real(vector(kk))*sqsmat(kk,jj)
         enddo
!         write(666,*) ii,jj,foruri
       enddo
       enddo

       return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
