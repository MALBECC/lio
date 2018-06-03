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
    real*8       :: Rho(2,2), S(2,2), outVec(2), criteria
    character*20 :: testResult
    integer      :: atomOrb(2)
    logical      :: testCond

    write(*,*) '- Mulliken Population Tests -'
    atomOrb(1) = 1
    atomOrb(2) = 2
    criteria   = 0.000000001d0

    ! Tests for null sqs matrix.
    Rho    = 1.0d0
    S      = 0.0d0
    outVec = 0.0d0
    testResult = "FAILED"

    call mulliken_calc( 2, 2, Rho, S, atomOrb, outVec)
    testCond = (abs(outVec(1)) < criteria).and.(abs(outVec(2)) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Obtain null vector from null S.'

    ! Tests for null rho matrix.
    Rho    = 0.0d0
    S      = 1.0d0
    outVec = 0.0d0
    testResult = "FAILED"

    call mulliken_calc(2 , 2, Rho, S, atomOrb, outVec)
    testCond = (abs(outVec(1)) < criteria).and.(abs(outVec(2)) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Obtain null vector from null RHO.'

    ! Tests for diagonal sqs matrix.
    Rho(1,1) = 1.0d0 ; Rho(1,2) = 0.0d0 ;
    Rho(2,1) = 0.0d0 ; Rho(2,2) = 1.0d0 ;

    S(1,1)   = 2.0d0 ; S(1,2)   = 2.0d0 ;
    S(2,1)   = 1.0d0 ; S(2,2)   = 1.0d0 ;

    outVec     = 0.0d0
    testResult = "FAILED"

    call mulliken_calc(2, 2, Rho, S, atomOrb, outVec)
    testCond = (abs(outVec(1)+2) < criteria).and.(abs(outVec(2)+1) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Usage of diagonal S matrix.'

    ! Tests for a sqs matrix with a negative value.
    Rho(1,1) = 2.0d0 ; Rho(1,2) = 1.0d0 ;
    Rho(2,1) =-1.0d0 ; Rho(2,2) = 3.0d0 ;

    S(1,1)   = 2.0d0 ; S(1,2)   = 2.0d0 ;
    S(2,1)   = 1.0d0 ; S(2,2)   = 1.0d0 ;

    outVec     = 0.0d0
    testResult = "FAILED"

    call mulliken_calc(2, 2, Rho, S, atomOrb, outVec)
    testCond = (abs(outVec(1)+6) < criteria).and.(abs(outVec(2)+2) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Usage of a negative element in S matrix.'

    ! Same test as before, but changing the value of M.
    outVec     = 0.0d0
    testResult = "FAILED"

    call mulliken_calc(1, 2, Rho, S,atomOrb, outVec)
    testCond = (abs(outVec(1)+6) < criteria).and.(abs(outVec(2)+2) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Usage of a single orbital.'

    return
end subroutine test_mulliken

subroutine test_lowdin()
    !this test should be redone and smat should be fed to lowdin charges
    implicit none
    real*8       :: Rho(2,2),SQS(2,2), outVec(2), criteria
    character*20 :: testResult
    integer      :: atomOrb(2)
    logical      :: testCond

    write(*,*) '- Löwdin Population Tests -'
    atomOrb(1) = 1
    atomOrb(2) = 2
    criteria   = 0.000000001d0
 
    ! Tests for null sqs matrix.
    Rho    = 1.0d0
    SQS    = 0.0d0
    outVec = 0.0d0
    testResult = "FAILED"

    call lowdin_calc( 2, 2, Rho, SQS, atomOrb, outVec)
    testCond = (abs(outVec(1)) < criteria).and.(abs(outVec(2)) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Obtain null vector from null SQS.'

    ! Tests for null rho matrix.
    Rho     = 0.0d0
    SQS     = 1.0d0
    outVec     = 0.0d0
    testResult = "FAILED"
  
    call lowdin_calc(2 , 2, Rho, SQS, atomOrb, outVec)
    testCond = (abs(outVec(1)) < criteria).and.(abs(outVec(2)) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Obtain null vector from null RHO.'

    ! Tests for diagonal sqs matrix.
    Rho(1,1) = 1.0d0 ; Rho(1,2) = 0.0d0 ;
    Rho(2,1) = 0.0d0 ; Rho(2,2) = 1.0d0 ;

    SQS(1,1) = 2.0d0 ; SQS(1,2) = 2.0d0 ;
    SQS(2,1) = 1.0d0 ; SQS(2,2) = 1.0d0 ;

    outVec     = 0.0d0
    testResult = "FAILED"

    call lowdin_calc(2, 2, Rho, SQS, atomOrb, outVec)
    testCond = (abs(outVec(1)+6) < criteria).and.(abs(outVec(2)+3) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Usage of diagonal SQS matrix.'

    ! Tests for a sqs matrix with a negative value.
    Rho(1,1) = 2.0d0 ; Rho(1,2) = 1.0d0 ;
    Rho(2,1) =-1.0d0 ; Rho(2,2) = 3.0d0 ;

    SQS(1,1) = 2.0d0 ; SQS(1,2) = 2.0d0 ;
    SQS(2,1) = 1.0d0 ; SQS(2,2) = 1.0d0 ;

    outVec(:)  = 0.0d0
    testResult = "FAILED"
 
    call lowdin_calc(2, 2, Rho, SQS, atomOrb, outVec)
    testCond = (abs(outVec(1)+12) < criteria).and.(abs(outVec(2)+6) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Usage of a negative element in SQS matrix.'

    ! Same test as before, but changing the value of M.
    outVec(:)  = 0.0d0
    testResult = "FAILED"

    call lowdin_calc(2, 1, Rho, SQS, atomOrb, outVec)
    testCond = (abs(outVec(1)+8) < criteria).and.(abs(outVec(2)) < criteria) 
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Usage of a single orbital.'

    return 
end subroutine test_lowdin
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% Orbital Energy related Functions.                                  [OEF] %%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine test_degeneration()
    implicit none
    integer              :: M, nDeg, nOrb
    integer, allocatable :: nDegMO(:)
    real*8 , allocatable :: energies(:)
    character*20         :: testResult
    logical              :: testCond
 
    write(*,*) '- Get Degeneration Tests -'
   
    ! Test for a single orbital.
    M = 1
    allocate(energies(M), nDegMO(M))
    nOrb       = 1
    energies   = 1.0d0
    nDeg       = 0
    nDegMO     = 0
    testResult = "FAILED"

    call get_degeneration(energies, nOrb, M, nDeg, nDegMO)   
    testCond = (nDeg.eq.1).and.(nDegMO(1).eq.1)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Usage of a single orbital.'
    deallocate(energies, nDegMO)

    ! Test for similar but not equal energies.
    M = 10
    allocate(energies(M), nDegMO(M))
    energies    = 1.0d0 ; energies(1) = 2.0d0
    energies(2) = 2.0d0 ; energies(3) = 2.000011d0
    nDeg   = 0
    nDegMO = 0
    nOrb   = 3
    testResult = "FAILED"

    call get_degeneration(energies, nOrb, M, nDeg, nDegMO)
    testCond = (nDeg.eq.1).and.(nDegMO(1).eq.3)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Similar but non-equal energies.'

    ! Test for 2 degenerate orbitals.
    nOrb = 1
    testResult = "FAILED"

    call get_degeneration(energies, nOrb, M, nDeg, nDegMO)
    testCond = (nDeg.eq.2).and.(nDegMO(1).eq.2).and.(nDegMO(2).eq.1)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Two degenerated orbitals.'

    ! Test for more than 2 degenerate orbitals.
    nOrb = 4
    testResult = "FAILED"

    call get_degeneration(energies, nOrb, M, nDeg, nDegMO)
    testCond = (nDeg.eq.7).and.(nDegMO(1).eq.4).and.(nDegMO(2).eq.5).and.      &
               (nDegMO(3).eq.6).and.(nDegMO(4).eq.7).and.(nDegMO(5).eq.8).and. &
               (nDegMO(6).eq.9).and.(nDegMO(7).eq.10)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - More than two degenerated orbitals.'
    deallocate(energies, nDegMO)
end subroutine test_degeneration

subroutine test_softness()
    implicit none
    real*8 :: enAH, enAL, enBH, enBL, soft, criteria
    character*20 :: testResult
    logical      :: testCond

    write(*,*) '- Get Softness -'

    enAH = -1.0d0 ; enBH = -2.0d0
    enAL =  3.0d0 ; enBL =  4.0d0
    soft =  0.0d0 ; criteria = 0.000000001d0
    testResult = "FAILED"

    call get_softness(enAH, enAL, enBH, enAL, soft)
    testCond = ((soft-0.4d0) < criteria)
    if (testCond) testResult = "PASSED"
    write(*,*) testResult, ' - Softness properly calculated.'
end subroutine test_softness
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% Fukui function tests.                                              [FUK] %%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

subroutine test_fukui()
! TO-DO: Proper Fukui function test.
end subroutine test_fukui

subroutine test_fukuiOS()
! TO-DO: Proper open-shell Fukui function test.
end subroutine test_fukuiOS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
