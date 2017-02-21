!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% FOCK_COMMUTS.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains calc_fock_commuts, which calculates F' and [F',P'] for    !
! DIIS matrices fockm and FP_PFm.                                              !
! Input: F (fock), P (rho), X, Y                                               !
! F' = X^T * F * X   |   P' = Y^T * P * Y   | X = (Y^-1)^T                     !
! => [F',P'] = X^T * F * P * Y - Y^T * P * F * X = A - A^T                     !
! Where A = X^T * F * P * Y                                                    !
! Output: A (scratch), A^T (scratch1), F' (fock)                               !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine calc_fock_commuts(fock,rho,X,Y,scratch,scratch1,M)

    implicit none
    integer, intent(in)    :: M
    REAL*8,  intent(in)    :: rho(M,M),X(M,M),Y(M,M)
    REAL*8,  intent(inout) :: fock(M,M)
    REAL*8,  intent(out)   :: scratch(M,M),scratch1(M,M)
    integer :: i, j, k

    ! X^T * F = scratch^T
    scratch = 0
    do i = 1, M
    do j = 1, M
    do k = 1, i ! X is upper triangular.
        scratch(j,i) = scratch(j,i) + X(k,i)*fock(k,j)
    enddo
    enddo
    enddo
    
    ! Do *X for fockm.
    fock=0
    do i = 1, M
    do j = 1, M
    do k = 1, j ! X is upper triangular
        fock(i,j) = fock(i,j) + scratch(k,i)*X(k,j) 
    enddo
    enddo
    enddo

    ! * P = scratch1^T
    scratch1=0
    do i = 1, M
    do j = 1, M
    do k = 1, M
        scratch1(j,i) = scratch1(j,i) + scratch(k,i)*rho(k,j)
    enddo
    enddo
    enddo

    ! * Y = scratch = scratch1^1
    scratch=0
    do i = 1, M
    do j = 1, M
        do k = j, M ! Y is lower triangular.
            scratch(i,j) = scratch(i,j) + scratch1(k,i)*Y(k,j)
        enddo
        ! The next k loop will only use values of scratch1(first index) > j, 
        ! so at this point we can overwrite anything behind that and just  
        ! reuse scratch1 as the transpose.
        scratch1(j,i) = scratch(i,j)
    enddo
    enddo
end subroutine calc_fock_commuts
