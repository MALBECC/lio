!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
include "../properties.f90"
program test_properties
 
    ! Electronic Population Analysis.     [EPA]
    call test_lowdin() 
    call test_mulliken()

    ! Orbital energy related functions.   [OEF]
    call test_degeneration()
    call test_softness()

    ! Fukui function tests.               [FUK]
    call test_fukui()
    call test_fukuiOS() 

end program test_properties




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% Electronic Population Analysis.                                    [EPA] %%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine test_mulliken()
    implicit none
    real*8       :: rhomat(2,2),sqsmat(2,2), outvec(2), criteria
    character*20 :: testResult
    integer      :: atomorb(2)
    logical      :: testCond

    write(*,*) '- Löwdin Population Tests -'
    atomorb(1) = 1
    atomorb(2) = 2
    criteria   = 0.000000001d0

    ! Tests for null sqs matrix.
    rhomat     = 1.0d0
    sqsmat     = 0.0d0
    outvec     = 0.0d0
    testResult = "FAILED"

    call lowdin_calc( 2, 2, rhomat, sqsmat, atomorb, outvec)
    testCond = (abs(outvec(1)) < criteria).and.(abs(outvec(2)) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Obtain null vector from null SQS.'

    ! Tests for null rho matrix.
    rhomat     = 0.0d0
    sqsmat     = 1.0d0
    outvec     = 0.0d0
    testResult = "FAILED"

    call lowdin_calc(2 , 2, rhomat, sqsmat, atomorb, outvec)
    testCond = (abs(outvec(1)) < criteria).and.(abs(outvec(2)) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Obtain null vector from null RHO.'

    ! Tests for diagonal sqs matrix.
    rhomat(1,1) = 1.0d0 ; rhomat(1,2) = 0.0d0 ;
    rhomat(2,1) = 0.0d0 ; rhomat(2,2) = 1.0d0 ;

    sqsmat(1,1) = 2.0d0 ; sqsmat(1,2) = 2.0d0 ;
    sqsmat(2,1) = 1.0d0 ; sqsmat(2,2) = 1.0d0 ;

    outvec     = 0.0d0
    testResult = "FAILED"

    call lowdin_calc(2,2,rhomat,sqsmat,atomorb,outvec)
    testCond = (abs(outvec(1)+6) < criteria).and.(abs(outvec(2)+3) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Usage of diagonal SQS matrix.'

    ! Tests for a sqs matrix with a negative value.
    rhomat(1,1) = 2.0d0 ; rhomat(1,2) = 1.0d0 ;
    rhomat(2,1) =-1.0d0 ; rhomat(2,2) = 3.0d0 ;

    sqsmat(1,1) = 2.0d0 ; sqsmat(1,2) = 2.0d0 ;
    sqsmat(2,1) = 1.0d0 ; sqsmat(2,2) = 1.0d0 ;

    outvec(:)  = 0.0d0
    testResult = "FAILED"

    call lowdin_calc(2,2,rhomat,sqsmat,atomorb,outvec)
    testCond = (abs(outvec(1)+12) < criteria).and.(abs(outvec(2)+6) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Usage of a negative element in SQS matrix.'

    ! Same test as before, but changing the value of M.
    outvec(:)  = 0.0d0
    testResult = "FAILED"

    call lowdin_calc(1, 2,rhomat,sqsmat,atomorb,outvec)
    testCond = (abs(outvec(1)+8) < criteria).and.(abs(outvec(2)) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Usage of a negative element in SQS matrix.'

    return
end subroutine test_mulliken

subroutine test_lowdin()
    implicit none
    real*8       :: rhomat(2,2),sqsmat(2,2), outvec(2), criteria
    character*20 :: testResult
    integer      :: atomorb(2)
    logical      :: testCond

    write(*,*) '- Löwdin Population Tests -'
    atomorb(1) = 1
    atomorb(2) = 2
    criteria   = 0.000000001d0
 
    ! Tests for null sqs matrix.
    rhomat     = 1.0d0
    sqsmat     = 0.0d0
    outvec     = 0.0d0
    testResult = "FAILED"

    call lowdin_calc( 2, 2, rhomat, sqsmat, atomorb, outvec)
    testCond = (abs(outvec(1)) < criteria).and.(abs(outvec(2)) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Obtain null vector from null SQS.'

    ! Tests for null rho matrix.
    rhomat     = 0.0d0
    sqsmat     = 1.0d0
    outvec     = 0.0d0
    testResult = "FAILED"
  
    call lowdin_calc(2 , 2, rhomat, sqsmat, atomorb, outvec)
    testCond = (abs(outvec(1)) < criteria).and.(abs(outvec(2)) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Obtain null vector from null RHO.'

    ! Tests for diagonal sqs matrix.
    rhomat(1,1) = 1.0d0 ; rhomat(1,2) = 0.0d0 ;
    rhomat(2,1) = 0.0d0 ; rhomat(2,2) = 1.0d0 ;

    sqsmat(1,1) = 2.0d0 ; sqsmat(1,2) = 2.0d0 ;
    sqsmat(2,1) = 1.0d0 ; sqsmat(2,2) = 1.0d0 ;

    outvec     = 0.0d0
    testResult = "FAILED"

    call lowdin_calc(2, 2, rhomat, sqsmat, atomorb, outvec)
    testCond = (abs(outvec(1)+6) < criteria).and.(abs(outvec(2)+3) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Usage of diagonal SQS matrix.'

    ! Tests for a sqs matrix with a negative value.
    rhomat(1,1) = 2.0d0 ; rhomat(1,2) = 1.0d0 ;
    rhomat(2,1) =-1.0d0 ; rhomat(2,2) = 3.0d0 ;

    sqsmat(1,1) = 2.0d0 ; sqsmat(1,2) = 2.0d0 ;
    sqsmat(2,1) = 1.0d0 ; sqsmat(2,2) = 1.0d0 ;

    outvec(:)  = 0.0d0
    testResult = "FAILED"
 
    call lowdin_calc(2, 2, rhomat, sqsmat, atomorb, outvec)
    testCond = (abs(outvec(1)+12) < criteria).and.(abs(outvec(2)+6) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Usage of a negative element in SQS matrix.'

    ! Same test as before, but changing the value of M.
    outvec(:)  = 0.0d0
    testResult = "FAILED"

    call lowdin_calc(1, 2, rhomat, sqsmat, atomorb, outvec)
    testCond = (abs(outvec(1)+8) < criteria).and.(abs(outvec(2)) < criteria) 
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Usage of a negative element in SQS matrix.'

    return 
end subroutine test_lowdin
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% Orbital Energy related Functions.                                  [OEF] %%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine test_degeneration()
end subroutine test_degeneration

subroutine test_softness()
end subroutine test_softness

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% Fukui function tests.                                              [FUK] %%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

subroutine test_fukui()
end subroutine test_fukui

subroutine test_fukuiOS()
end subroutine test_fukuiOS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
