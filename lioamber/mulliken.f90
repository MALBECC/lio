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
!         q(kk)=0. !real(q0(kk))
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
        write(UID,*)
        write(UID,300)
        write(UID,301)
        write(UID,302)
        write(UID,303)
        write(UID,304)
       do kk=1,N
         write(UID,305) kk,q0(kk),q(kk)
       enddo
       write(UID,306)
       write(UID,*)

       call g2g_timer_stop('mulliken_write')

 300   FORMAT(8x,"╔═════════════════&
       ════════════════╗")
 301   FORMAT(8x,"║   MULLIKEN POPULATION ANALYSIS  ║")
 302   FORMAT(8x,"╠════════╦═══════════╦════════════╣")
 303   FORMAT(8x,"║ ATOM # ║ ATOM TYPE ║ POPULATION ║")
 304   FORMAT(8x,"╠════════╬═══════════╬════════════╣")
 305   FORMAT(8x,"║",2x,i3,3x,"║"3x,i3,5x,"║",1x,F10.7,1x,"║")
 306   FORMAT(8x,"╚════════╩═══════════╩════════════╝")
       return
       end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

