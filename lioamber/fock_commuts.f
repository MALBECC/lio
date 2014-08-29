!-----------------------------------------------------------------------------------------
! Calculates F' and [F',P'] for diis matrices fockm and FP_PFm
!
! Input: F (fock), P (rho), X, Y
!
! F' = X^T * F * X
! P' = Y^T * P * Y
! X = (Y^-1)^T
! => [F',P'] = X^T * F * P * Y - Y^T * P * F * X
! =>         = A - A^T
! where A = X^T * F * P * Y
!
! Output: A (scratch2), A^T (scratch1), F' (fock)
!-----------------------------------------------------------------------------------------
      subroutine calc_fock_commuts(fock,rho,X,Y,scratch,
     >                             scratch1,scratch2,M)

          integer, intent(in)    :: M
          REAL*8,  intent(in)    :: rho(M,M),X(M,M),Y(M,M)
          REAL*8,  intent(inout) :: fock(M,M)
          REAL*8,  intent(out)   :: scratch(M,M),scratch1(M,M)
          REAL*8,  intent(out)   :: scratch2(M,M)

          integer i,j,k

!---------------------------------------------------------------------
! X^T * F = scratch^T
!---------------------------------------------------------------------
          scratch=0
          do i=1,M
            do j=1,M
              ! X is upper triangular
              do k=1,i
                scratch(j,i)=scratch(j,i)+X(k,i)*fock(k,j)
              enddo
            enddo
          enddo
!---------------------------------------------------------------------
! do * X for fockm
!---------------------------------------------------------------------
          fock=0
          do i=1,M
            do j=1,M
              ! X is upper triangular
              do k=1,j
                fock(i,j)=fock(i,j)+scratch(k,i)*X(k,j) 
              enddo
            enddo
          enddo
!---------------------------------------------------------------------
! * P = scratch1^T
!---------------------------------------------------------------------
          scratch1=0
          do i=1,M
            do j=1,M
              do k=1,M
                scratch1(j,i)=scratch1(j,i)+scratch(k,i)*rho(k,j)
              enddo
            enddo
          enddo
!---------------------------------------------------------------------
! * Y = scratch2 = scratch1^1
!---------------------------------------------------------------------
          scratch2=0
          do i=1,M
            do j=1,M
              ! Y is lower triangular
              do k=j,M
                scratch2(i,j)=scratch2(i,j)+scratch1(k,i)*Y(k,j)
              enddo
              scratch1(j,i)=scratch2(i,j)
            enddo
         enddo
      end
