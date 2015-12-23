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
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: UID,N
        INTEGER, INTENT(IN) :: q0(N)
        REAL*8, INTENT(IN)  :: q(N)
        REAL*8              :: qtotal
        INTEGER             :: kk
        CALL g2g_timer_start('mulliken_write')
!
        qtotal=0.d0
        WRITE(UID,*)
        WRITE(UID,300)
        WRITE(UID,301)
        WRITE(UID,302)
        WRITE(UID,303)
        WRITE(UID,304)
        DO kk=1,N
         qtotal=qtotal+q(kk)
         WRITE(UID,305) kk,q0(kk),q(kk)
        ENDDO
        WRITE(UID,306)
        WRITE(UID,307) qtotal
        WRITE(UID,308)
        WRITE(UID,*)

        CALL g2g_timer_stop('mulliken_write')

 300   FORMAT(8x,"╔═════════════════&
       ════════════════╗")
 301   FORMAT(8x,"║   MULLIKEN POPULATION ANALYSIS  ║")
 302   FORMAT(8x,"╠════════╦═══════════╦════════════╣")
 303   FORMAT(8x,"║ ATOM # ║ ATOM TYPE ║ POPULATION ║")
 304   FORMAT(8x,"╠════════╬═══════════╬════════════╣")
 305   FORMAT(8x,"║",2x,i3,3x,"║"3x,i3,5x,"║",1x,F10.7,1x,"║")
 306   FORMAT(8x,"╚════════╬═══════════╬════════════╣")
 307   FORMAT(8x,"         ║   TOTAL   ║",1x,F10.7,1x,"║") 
 308   FORMAT(8x,"         ╚═══════════╩════════════╝")
       return
       END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

