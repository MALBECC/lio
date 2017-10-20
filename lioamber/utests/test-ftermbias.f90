!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
include "../fterm_biaspot.f"
program test_ftermbias
    implicit none
    real*8       :: sqsmat(3,3),outmat(3,3), weight, criteria
    integer      :: vector(3)
    logical      :: testCond
    character*20 :: testResult

    criteria = 0.000000001d0
    write(*,*) '- FTerm Bias Potential Tests -'

    ! Tests for null SQS matrix.
    sqsmat = 0.0d0
    vector = 1
    weight = 1.0d0
    outmat = 0.0d0
    testResult = "FAILED"
    
    call fterm_biaspot(3, sqsmat, vector, weight, outmat)
    testCond = (maxval(abs(outmat)) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Obtain null matrix from null SQS.'
  
    ! Tests for null vector.
    sqsmat(1,1) = 1.0d0 ; sqsmat(1,2) = 3.0d0 ; sqsmat(1,3) = 2.0d0
    sqsmat(2,1) = 5.0d0 ; sqsmat(2,2) = 6.0d0 ; sqsmat(2,3) = 3.0d0
    sqsmat(3,1) = 4.0d0 ; sqsmat(3,2) = 2.0d0 ; sqsmat(3,3) = 1.0d0

    vector = 0
    outmat = 0.0d0
    testResult = "FAILED"

    call fterm_biaspot(3, sqsmat, vector, weight, outmat)
    testCond = (maxval(abs(outmat)) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Obtain null matrix from null vector.'

    ! Tests for null weight.
    vector = 1
    weight = 0.0d0
    outmat = 0.0d0
    testResult = "FAILED"

    call fterm_biaspot(3, sqsmat, vector, weight, outmat)
    testCond = (maxval(abs(outmat)) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Obtain null matrix from null weight.'

    ! Tests for diagonal SQS matrix.
    sqsmat      = 0.0d0
    sqsmat(1,1) = 2.0d0
    sqsmat(2,2) = 3.0d0
    sqsmat(3,3) = 4.0d0

    weight = 1.0d0
    outmat = 0.0d0

    call fterm_biaspot(3, sqsmat, vector, weight, outmat)
    testCond = (abs(outmat(1,1) - 4 ) < criteria).and. &
               (abs(outmat(2,2) - 9 ) < criteria).and. &
               (abs(outmat(3,3) - 16) < criteria).and. &
               (minval(abs(outmat)) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Test for diagonalized SQS matrix.'

    ! Tests diagonal SQS with a 0 in vector.
    vector    = 1
    vector(3) = 0
    weight    = 2.0d0
    outmat    = 0.0d0
   
    call fterm_biaspot(3, sqsmat, vector, weight, outmat)
    testCond = (abs(outmat(1,1) - 8 ) < criteria).and. &
               (abs(outmat(2,2) - 18) < criteria).and. &
               (minval(abs(outmat)) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Test for (1, 1, 0) vector.'

    ! Tests for all 1.
    sqsmat = 1.0d0
    vector = 1
    weight = 1.0d0

    call fterm_biaspot(3, sqsmat, vector, weight, outmat)
        testCond = (abs(outmat(1,1) - 11) < criteria).and. &
                   (abs(outmat(2,2) - 21) < criteria).and. &
                   (minval(abs(outmat-3)) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Test for all 1 SQS, vector, and weight.'

    return
end program test_ftermbias
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
