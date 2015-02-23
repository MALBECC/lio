!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! PACKED STORAGE
!------------------------------------------------------------------------------!
!
! This two subroutines are intended to facilitate the passage
! from a matrix and a vector with the same information (the
! vector has it in what is known as "packed storage" or SP).
! The matrix must be symetric for packed storage to be used,
! and this routines only support real values.
!
! NOTE: If Vector is originally of size > NM*(NM+1)/2, then
! calling the subroutine with "Vector(i0)" will make it
! work with the positions of the original vector that go
! from i0 to i0+NM*(NM+1)/2-1.
! (something similar may work with the variable 'Matrix'
! but it will probably be more intrincated)
!
!-----------------------!
! STORAGE SCHEME UPLO=U
!-----------------------!
!  | 1  7  8 |
!  |(7) 2  9 |   <=>   ( 1 , 7 , 2 , 8 , 9 , 3 )
!  |(8)(9) 3 |
!
!-----------------------!
! STORAGE SCHEME UPLO=L
!-----------------------!
!  | 1 (7)(8)|
!  | 7  2 (9)|   <=>   ( 1 , 7 , 8 , 2 , 9 , 3 )
!  | 8  9  3 |
!
!----------------------------------------------------------!
! For more information, visit:
!   (*) http://www.netlib.org/lapack/lug/node123.html
!   (*) http://www.netlib.org/lapack/lug/node24.html
!
! 04/2014 || F.F.R
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE sprepack(UPLO,NM,Vector,Matrix)
       IMPLICIT NONE
       CHARACTER(LEN=1)   :: UPLO
       INTEGER,INTENT(IN) :: NM
       REAL*8,INTENT(OUT) :: Vector(NM*(NM+1)/2)
       REAL*8,INTENT(IN)  :: Matrix(NM,NM)
       INTEGER            :: ii,jj,idx
!
!------------------------------------------------------------------------------!
!
       IF (UPLO.EQ.'U') THEN
         DO jj=1,NM
         DO ii=1,jj
           idx=ii+(jj*(jj-1)/2)
           Vector(idx)=Matrix(ii,jj)
         ENDDO
         ENDDO

       ELSE IF (UPLO.EQ.'L') THEN
         DO ii=1,NM
         DO jj=1,ii
           idx=ii+(2*NM-jj)*(jj-1)/2
           Vector(idx)=Matrix(ii,jj)
         ENDDO
         ENDDO

       ELSE
         PRINT*,'NOT GOOD INPUT FOR UPLO'

       ENDIF
!
!------------------------------------------------------------------------------!
       RETURN;END SUBROUTINE


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
