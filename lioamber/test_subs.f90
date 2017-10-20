!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% TEST_SUBS.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains subroutines intended to be used in debug.                 !
! Subroutines included are:                                                    !
! * seek_nan (Seeks NAN in vectors)                                            !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine seek_nan(vecToTest, vecStart, vecEnd, phrase)
    implicit none
    real*8           , intent(in) :: vecToTest(*)     ! Vector to analize.
    integer          , intent(in) :: vecStart, vecEnd ! Vector range to analize.
    character (len=*), intent(in) :: phrase           ! Output phrase for NaN.
    integer :: iNick

    if (vecStart .gt. vecEnd) then
        write(*,*) "Error: vector start index greater than end index."
        write(*,*) phrase
        stop
    endif

    do iNick = vecStart, vecEnd
        if (vecToTest(iNick) .ne. vecToTest(iNick)) then
            write(*,*) "NaN found in: ", phrase, iNick
            stop
        end if
    enddo
endsubroutine seek_nan
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
