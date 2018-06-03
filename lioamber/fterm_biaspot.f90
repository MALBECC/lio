!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% FTERM_BIASPOT.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! TO-DO: EXPLANATION.
! This subroutine is used by the ehrenfest module only
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fterm_biaspot(M,sqsmat,vector,weight,outmat)
    implicit none
    integer, intent(in)    :: M, vector(M)
    real*8 , intent(in)    :: weight, sqsmat(M,M)
    real*8 , intent(inout) :: outmat(M,M)

    real*8  :: newterm, foruri
    integer :: ii, jj, kk

    do ii = 1, M
    do jj = 1, M
    do kk = 1, M
        newterm       = sqsmat(ii,kk) * weight * sqsmat(kk,jj)
        outmat(ii,jj) = outmat(ii,jj) + newterm*real(vector(kk))
    enddo
    enddo
    enddo

    return
end subroutine fterm_biaspot
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
