!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% PACKED_STORAGE.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This two subroutines are intended to facilitate the passage from a matrix and!
! a vector with the same information (the vector has it in what is known as    !
! "packed storage" or SP). The matrix must be symetric for packed storage to   !
! be used, and this routines only support real values.                         !
! Subroutines included are:                                                    !
! * spunpack                                                                   !
! * spunpack_rho                                                               !
! * sprepack                                                                   !
! * spunpack_rtc                                                               !
! * sprepack_ctr                                                               !
!                                                                              !
! NOTE: If Vector is originally of size > NM*(NM+1)/2, then calling the        !
! subroutine with "Vector(i0)" will make it work with the positions of the     !
! original vector that go from i0 to i0+NM*(NM+1)/2-1. (something similar may  !
! work with the variable 'Matrix' but it will probably be more intrincated).   !
!                                                                              !
! STORAGE SCHEME UPLO=U                                                        !
!  | 1  7  8 |                                                                 !
!  |(7) 2  9 |   <=>   ( 1 , 7 , 2 , 8 , 9 , 3 )                               !
!  |(8)(9) 3 |                                                                 !
! STORAGE SCHEME UPLO=L                                                        !
!  | 1 (7)(8)|                                                                 !
!  | 7  2 (9)|   <=>   ( 1 , 7 , 8 , 2 , 9 , 3 )                               !
!  | 8  9  3 |                                                                 !
!                                                                              !
! For more information, visit:                                                 !
!   (*) http://www.netlib.org/lapack/lug/node123.html                          !
!   (*) http://www.netlib.org/lapack/lug/node24.html                           !
! 04/2014 || F.F.R                                                             !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine spunpack(UPLO, NM, Vector, Matrix)
    implicit none
    character(len=1)    :: UPLO
    integer,intent(in)  :: NM
    real*8 ,intent(in)  :: Vector(NM*(NM+1)/2)
    real*8 ,intent(out) :: Matrix(NM,NM)
    integer             :: ii, jj, idx

    if (UPLO.eq.'U') then
        do jj = 1, NM
        do ii = 1, jj
            idx = ii + (jj*(jj-1)/2)
            Matrix(ii,jj) = Vector(idx)
            Matrix(jj,ii) = Vector(idx)
        enddo
        enddo

    else if (UPLO.eq.'L') then
        do ii = 1, NM
        do jj = 1, ii
            idx = ii + (2*NM-jj)*(jj-1)/2
            Matrix(ii,jj) = Vector(idx)
            Matrix(jj,ii) = Vector(idx)
        enddo
        enddo

    else
        write(*,*) 'spunpack: Wrong value in UPLO.'
    endif

    return
end subroutine spunpack

subroutine spunpack_rho(UPLO, NM, Vector, Matrix)
    implicit none
    character(len=1)     :: UPLO
    integer, intent(in)  :: NM
    real*8 , intent(in)  :: Vector(NM*(NM+1)/2)
    real*8 , intent(out) :: Matrix(NM,NM)
    integer              :: ii, jj, idx

    if (UPLO.eq.'U') then
        do jj = 1, NM
            do ii = 1, jj - 1
               idx = ii + (jj*(jj-1)/2)
               Matrix(ii,jj) = Vector(idx)/2
               Matrix(jj,ii) = Vector(idx)/2
            enddo
               idx = jj + (jj*(jj-1)/2)
               Matrix(jj,jj) = Vector(idx)
        enddo
    else if (UPLO.eq.'L') then
        do ii = 1, NM
            do jj = 1, ii - 1
                idx = ii+(2*NM-jj)*(jj-1)/2
                Matrix(ii,jj) = Vector(idx)/2
                Matrix(jj,ii) = Vector(idx)/2
            enddo
            idx = jj + (2*NM-jj)*(jj-1)/2
            Matrix(jj,jj) = Vector(idx)
        enddo
    else
        write(*,*) 'spunpack_rho: Wrong value in UPLO.'
    endif

    return
end subroutine spunpack_rho

subroutine sprepack(UPLO, NM, Vector, Matrix)
    implicit none
    character(len=1)    :: UPLO
    integer,intent(in)  :: NM
    real*8 ,intent(in)  :: Matrix(NM, NM)
    real*8 ,intent(out) :: Vector(NM*(NM+1)/2)
    integer             :: ii, jj, idx

    if (UPLO.eq.'U') then
        do jj = 1, NM
        do ii = 1, jj
            idx = ii + (jj*(jj-1)/2)
            Vector(idx) = Matrix(ii,jj)
        enddo
        enddo
    else if (UPLO.eq.'L') then
        do ii = 1, NM
        do jj = 1, ii
           idx = ii + (2*NM-jj)*(jj-1)/2
           Vector(idx) = Matrix(ii,jj)
        enddo
        enddo
    else
         write(*,*) 'sprepack: Wrong value in UPLO.'
    endif

    return
end subroutine sprepack

subroutine spunpack_rtc(UPLO, NM, Vector, Matrix)
    implicit none
    character(len=1)        :: UPLO
    integer   , intent(in)  :: NM
    real*8    , intent(in)  :: Vector(NM*(NM+1)/2)
#ifdef TD_SIMPLE
    complex*8 , intent(out) :: Matrix(NM,NM)
#else
    complex*16, intent(out) :: Matrix(NM,NM)
#endif
    integer                 :: ii, jj, idx

    if (UPLO.eq.'U') then
        do jj = 1, NM
            do ii = 1, jj - 1
                idx = ii + (jj*(jj-1)/2)
                Matrix(ii,jj) = cmplx(Vector(idx),0.0D0)
                Matrix(ii,jj) = Matrix(ii,jj)*0.50D0
                Matrix(jj,ii) = Matrix(ii,jj)
            enddo
            idx = ii + (ii*(ii-1)/2)
            Matrix(ii,ii) = cmplx(Vector(idx),0.0D0)
        enddo
    else if (UPLO.eq.'L') then
        do ii=1,NM
            do jj = 1, ii - 1
                idx = ii + (2*NM-jj)*(jj-1)/2
                Matrix(ii,jj) = cmplx(Vector(idx),0.0D0)
                Matrix(ii,jj) = Matrix(ii,jj)*0.50D0
            enddo
            Matrix(ii,ii) = Vector(ii+(2*NM-ii)*(ii-1)/2)
            do jj = ii + 1, NM
                idx = jj + (2*NM-ii)*(ii-1)/2
                Matrix(ii,jj) = cmplx(Vector(idx),0.0D0)
                Matrix(ii,jj) = Matrix(ii,jj)*0.50D0
            enddo
        enddo
    else
        write(*,*) 'spunpack_rtc: Wrong value in UPLO.'
    endif

    return
end subroutine spunpack_rtc

subroutine sprepack_ctr(UPLO,NM,Vector,Matrix)
    implicit none
    character(len=1)        :: UPLO
    integer   ,intent(in)   :: NM
    real*8    , intent(out) :: Vector(NM*(NM+1)/2)
#ifdef TD_SIMPLE
    complex*8 , intent(in)  :: Matrix(NM,NM)
#else
    complex*16, intent(in)  :: Matrix(NM,NM)
#endif
    integer                 :: ii, jj, idx

    if (UPLO.eq.'U') then
        do jj = 1, NM
        do ii = 1, jj
           idx = ii + (jj*(jj-1)/2)
           if(ii.eq.jj) then
              Vector(idx) = real(Matrix(ii,jj))
           else
              Vector(idx) = real(Matrix(ii,jj))
              Vector(idx) = Vector(idx)*2.0D0
           endif
        enddo
        enddo
    else if (UPLO.eq.'L') then
        do ii = 1, NM
        do jj = 1, ii
            idx = ii + (2*NM-jj)*(jj-1)/2
            if(ii.eq.jj) then
                Vector(idx) = real(Matrix(ii,jj))
            else
                Vector(idx) = real(Matrix(ii,jj))
                Vector(idx) = Vector(idx)*2.0D0
           endif
        enddo
        enddo
    else
        write(*,*) 'spunpack_rtc: Wrong value in UPLO.'
    endif

    return
end subroutine sprepack_ctr
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
