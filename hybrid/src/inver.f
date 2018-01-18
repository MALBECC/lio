
      SUBROUTINE INVER(A,B,N,NDIM,INFO)

      IMPLICIT NONE
      INTEGER N,NDIM,INFO
      DOUBLE PRECISION A(NDIM,NDIM),B(NDIM,NDIM),X

C Routine to obtain the inverse of a general, nonsymmetric
C square matrix, by using the Numerical Recipes algorithm.                
C The matrix A is of order N, but is defined with a 
C size NDIM >= N.
C
C P. Ordejon, June 2003
C
C **** INPUT *****
C A(NDIM,NDIM)   real*8      Matrix to be inverted
C N              integer     Size of the matrix
C NDIM           integer     Defined size of matrix
C **** OUTPUT ****
C B(NDIM,NDIM)   real*8      Inverse of A
C INFO           integer     If inversion was unsucessful, 
C                            INFO <> 0 is returned
C                            Is successfull, INFO = 0 returned
C ****************

      DOUBLE PRECISION WORK(N),C,ERROR,DELTA,TOL
      INTEGER IPIV(N)
      INTEGER I,J,K

      TOL=1.0D-4
      INFO = 0

C Try simple, old algorithm:

      DO I=1,N
        DO J=1,N
          B(I,J)=A(I,J)
        ENDDO
      ENDDO
      DO I=1,N
        IF (B(I,I) .EQ. 0.0D0) THEN
          INFO = 1
          RETURN
        ENDIF
        X=B(I,I)
        B(I,I)=1.0d0
        DO J=1,N
          B(J,I)=B(J,I)/X
        ENDDO
        DO K=1,N
          IF ((K-I).NE.0) THEN 
            X=B(I,K)
            B(I,K)=0.0d0
            DO J=1,N
              B(J,K)=B(J,K)-B(J,I)*X
            ENDDO
          ENDIF
        ENDDO
      ENDDO

C CHECK THAT THE INVERSE WAS CORRECTLY CALCULATED

      ERROR=0.0D0
      DO I=1,N
      DO J=1,N
        C=0.0D0
        DO K=1,N
          C=C+A(I,K)*B(K,J)
        ENDDO
        DELTA=0.0D0
        IF (I.EQ.J)  DELTA=1.0D0
        ERROR=ERROR+DABS(C-DELTA)
        C=0.0D0
        DO K=1,N
          C=C+B(I,K)*A(K,J)
        ENDDO
        DELTA=0.0D0
        IF (I.EQ.J)  DELTA=1.0D0
        ERROR=ERROR+DABS(C-DELTA)
      ENDDO
      ENDDO

      IF (ERROR/N .GT. TOL) THEN
        WRITE(6,'(a,f10.5)')
     .    'inver:  INVER unsuccessful, ERROR = ', ERROR
        INFO = 1
        RETURN
      ENDIF

      INFO = 0
      RETURN
      END

