!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!       pure subroutine lowdinpop(M,N,rhomat,sqsmat,atomorb,atomicq)
       subroutine lowdinpop(M,N,rhomat,sqsmat,atomorb,atomicq)
       implicit none
       integer,intent(in)   :: M,N
       real*8,intent(in)    :: rhomat(M,M)
       real*8,intent(in)    :: sqsmat(M,M)
       integer,intent(in)   :: atomorb(M)
       real*8,intent(inout) :: atomicq(N)

       real*8  :: newterm
       integer :: natom
       integer :: ii,jj,kk
!------------------------------------------------------------------------------!

       do kk=1,M
         natom=atomorb(kk)
         do ii=1,M
         do jj=1,M
           newterm=sqsmat(kk,ii)*rhomat(ii,jj)*sqsmat(jj,kk)
           atomicq(natom)=atomicq(natom)-newterm
         enddo
         enddo
       enddo

!------------------------------------------------------------------------------!
       return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
