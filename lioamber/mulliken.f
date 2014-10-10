!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine mulliken_calc(N,M,RealRho,Smat,NofM,q0,q)
!------------------------------------------------------------------------------!
!
! Rho must be written in atomic orbitals.
!
!------------------------------------------------------------------------------!
       implicit none
       integer,intent(in)   :: N,M,NofM(M)
       real*8,intent(in)    :: RealRho(M,M),Smat(M,M)
       integer,intent(in)   :: q0(N)
       real*8,intent(out)   :: q(N)
       integer              :: ii,jj,kk
       real*8               :: qe
       call g2g_timer_start('mulliken')
!
       do kk=1,N
         q(kk)=real(q0(kk))
       enddo

       do ii=1,M
       do jj=1,M
         qe=RealRho(ii,jj)*Smat(ii,jj)
         q(NofM(ii))=q(NofM(ii))-qe
       enddo
       enddo
!
       call g2g_timer_stop('mulliken')
       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE mulliken_write(UID,N,q0,q)
!------------------------------------------------------------------------------!
       implicit none
       integer,intent(in) :: UID,N
       integer,intent(in) :: q0(N)
       real*8,intent(in)  :: q(N)
       integer            :: kk
       call g2g_timer_start('mulliken_write')
!
       write(UID,300)
       write(UID,200) 'MULLIKEN POPULATION ANALYSIS'
       write(UID,200)
       write(UID,201) 'ATOM #','ATOM TYPE','POPULATION'
       do kk=1,N
         write(UID,202) kk,q0(kk),q(kk)
       enddo
       write(UID,200)
!
       call g2g_timer_stop('mulliken_write')
 200   format(A)
 201   format(A,4x,A,4x,A)
 202   format(I3,9X,I3,6X,F14.8)
 300   format('##################################################'
     > ,'##########')
       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

