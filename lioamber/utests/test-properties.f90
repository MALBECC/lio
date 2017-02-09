!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
program test_properties

    call test_lowdin()

end program test_properties

subroutine test_lowdin()
    implicit none
    real*8   :: rhomat(2,2),sqsmat(2,2), outvec(2), criteria
    integer  :: atomorb(2)

    write(*,*) ''
    write(*,*) ''
    atomorb(1) = 1
    atomorb(2) = 2
    criteria   = 0.000000001d0
 
    ! Tests for null sqs matrix.
    rhomat     = 1.0d0
    sqsmat     = 0.0d0
    outvec     = 0.0d0

    call lowdin_calc( 2, 2, rhomat, sqsmat, atomorb, outvec)
    write(*,*) 'Obtain null vector from null sqs:'
    write(*,*) outvec

    if ((abs(outvec(1)) < criteria).and.(abs(outvec(2)) < criteria)) & 
    write(*,*) "PASSED"

    ! Tests for null rho matrix.
    rhomat = 0.0d0
    sqsmat = 1.0d0
    outvec = 0.0d0
  
    call lowdin_calc(2,2,rhomat,sqsmat,atomorb,outvec)
    write(*,*) 'Obtain null vector from null rho:'
    write(*,*) outvec
    if ((abs(outvec(1)) < criteria).and.(abs(outvec(2)) < criteria)) &
    write(*,*) "PASSED"

    ! Tests for diagonal rho matrix.
    rhomat(1,1) = 1.0d0 ; rhomat(1,2) = 0.0d0 ;
    rhomat(2,1) = 0.0d0 ; rhomat(2,2) = 1.0d0 ;

    sqsmat(1,1) = 2.0d0 ; sqsmat(1,2) = 2.0d0 ;
    sqsmat(2,1) = 1.0d0 ; sqsmat(2,2) = 1.0d0 ;

    outvec = 0.0d0

    call lowdin_calc(2,2,rhomat,sqsmat,atomorb,outvec)
    write(*,*) 'expected result is (-6, -3) :'
    write(*,*) outvec
    if ((abs(outvec(1)+6) < criteria).and.(abs(outvec(2)+3) < criteria)) &
    write(*,*) "PASSED"

    ! Tests for a rho matrix with a negative value.
    rhomat(1,1) = 2.0d0 ; rhomat(1,2) = 1.0d0 ;
    rhomat(2,1) =-1.0d0 ; rhomat(2,2) = 3.0d0 ;

    sqsmat(1,1) = 2.0d0 ; sqsmat(1,2) = 2.0d0 ;
    sqsmat(2,1) = 1.0d0 ; sqsmat(2,2) = 1.0d0 ;

    outvec(:)=0.0d0
 
    call lowdin_calc(2,2,rhomat,sqsmat,atomorb,outvec)
    write(*,*) 'expected result is (-12, -6) :'
    write(*,*) outvec(1),outvec(2)
    if ((abs(outvec(1)+12) < criteria).and.(abs(outvec(2)+6) < criteria)) &
    write(*,*) "PASSED"

    ! Same test as before, but changing the value of M.
    outvec(:)=0.0d0

    call lowdin_calc(1, 2,rhomat,sqsmat,atomorb,outvec)
    write(*,*) 'expected result is (-8,0) :'
    write(*,*) outvec(1),outvec(2)
    if ((abs(outvec(1)+8) < criteria).and.(abs(outvec(2)) < criteria)) &
    write(*,*) "PASSED"

    return 
end subroutine test_lowdin
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
