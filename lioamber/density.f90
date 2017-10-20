!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DENSTITY.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains the routine density, which performs the construction of   !
! density matrix from the orbital coefficients.                                !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine density(M, NCO, X, C)
    use garcha_mod, only : OPEN   
    ! M:    number of atomic basis functions.
    ! NCO:  number of occupied orbitals.
    ! X:    matrix containing Fock coefficients.
    ! C:    matrix containing output density matrix.
    implicit none
    integer,intent(in)  :: M, NCO
    real*8, intent(in)  :: X(M, M)
    real*8, intent(out) :: C(M, M)

    real*8, allocatable :: A(:, :)
    real*8              :: factor
    integer             :: DOSM, i, j
#ifdef CUBLAS
    integer*8 devPtrA, devPtrC
    integer   sizeof_real
    parameter(sizeof_real=8)
    integer   stat
    external  CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
    external  CUBLAS_SHUTDOWN, CUBLAS_ALLOC,CUBLAS_FREE
    integer   CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
    integer   CUBLAS_INIT
#endif

    ! This factor take into account occupation for closed-shell
    ! systems.
    factor = 2.0D0 
    if (OPEN) factor = 1.0D0

    allocate(A(M,NCO))
    do i=1, M
        do j=1, NCO
             A(i, j) = X(i, j)
        enddo
    enddo
    C = 0.0D0

#ifdef CUBLAS
    ! Allocates matrix A on the device.
    stat = CUBLAS_ALLOC(M*NCO, sizeof_real, devPtrA)
    if (stat.NE.0) then
        write(*,*) "density: A-Matrix allocation failed."
        call CUBLAS_SHUTDOWN
        stop
    endif

    ! Allocate matrix C on the device
    stat = CUBLAS_ALLOC(M*M, sizeof_real, devPtrC)
    if (stat.NE.0) then
        write(*,*) "density: C-Matrix allocation failed."
        call CUBLAS_SHUTDOWN
        stop
    endif
 
    ! Copy A in the device
    stat = CUBLAS_SET_MATRIX(M, NCO, sizeof_real, A, M, devPtrA, M)
    if (stat.NE.0) then
        write(*,*) "matrix copy to the device failed -density"
        call CUBLAS_SHUTDOWN
        stop 
    endif

    ! Peforms the matrix multiplication, obtaining C as A*A.
    call CUBLAS_DGEMM('N', 'T', M, M, NCO, factor, devPtrA, M, devPtrA, M, &
                      0.0D0, devPtrC, M)

    ! Copies C to host.
    stat = CUBLAS_GET_MATRIX(M, M, sizeof_real, devPtrC, M, C, M)

    ! Finalizes variables.
    call CUBLAS_FREE(devPtrA)
    call CUBLAS_FREE(devPtrC)

#else
    ! Obtains C as A*A.
    call DGEMM('N', 'T', M, M, NCO, factor, A, M, A, M, 0.0D0, C, M)
#endif
    
    ! Multiplies by 2 all non diagonal elements.
    do i=1,M
        do j=1,i-1
            C(i,j) = 2.0D0 * C(i,j)
        enddo
        do j=i+1,M
            C(i,j) = 2.0D0 * C(i,j)
        enddo
    enddo

    deallocate(A)
    return
end subroutine density
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
