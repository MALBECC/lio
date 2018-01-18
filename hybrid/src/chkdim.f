      SUBROUTINE CHKDIM (SUB,VAR,ND,N,IOPT)
C **********************************************************************
C Checks a dimension parameter, warning and stopping if it is too small
C Written by J.M.Soler. 1987.
C ************ INPUT ***************************************************
C CHARACTER SUB*(*) : Name of calling subroutine
C CHARACTER VAR*(*) : Name of dimension parameter
C INTEGER   ND      : Actual value of dimension parameter
C INTEGER   N       : Required value of dimension parameter
C INTEGER   IOPT    : Option switch: IOPT=0 => Require that ND.EQ.N
C                                    IOPT=1 => Require that ND.GE.N
C **********************************************************************
C
C  Modules
C
      use precision
      use sys

      implicit none

      character(len=*), intent(in) ::  SUB ,VAR
      integer, intent(in)          ::  nd, n, iopt

      character(len=132) message

      IF ( IOPT.EQ.0 ) THEN
        IF ( ND.EQ.N ) RETURN
        WRITE (message,'(5A,I8,A,I8)') 'chkdim: ERROR: In ', SUB,
     .        ', dimension ',VAR,' =',ND,'. It must be exactly ',N
        call die(message)
      ELSE
        IF ( ND.GE.N ) RETURN
        WRITE (message,'(5A,I8,A,I8)') 'chkdim: ERROR: In ', SUB,
     .        ', dimension ',VAR,' =',ND,'. It must be at least ',N
        call die(message)
      ENDIF
      END



      SUBROUTINE CHKDIME (ND,N,OVERFLOW,NM)
      LOGICAL OVERFLOW
      integer nd, n, nm

      NM = MAX(N,NM)
      IF ( ND.GE.N ) RETURN
      OVERFLOW = .TRUE.
      RETURN
      END


