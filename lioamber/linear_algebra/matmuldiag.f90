!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% MATMULDIAG.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs C = A*B for square matrices.                                        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine matmuldiag(A, B, C, M)
    implicit none
    integer, intent(in)  :: M
    real*8,  intent(in)  :: A(M,M), B(M,M)
    real*8,  intent(out) :: C(M,M)
    integer :: i, k

    C = 0.0D0
    do k = 1, M
    do i = 1, M
        C(i,i) = C(i,i) + A(i,k)*B(k,i)
    enddo
    enddo

    return
end subroutine matmuldiag
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
