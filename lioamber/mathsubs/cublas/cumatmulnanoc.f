!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! REDUNDANT: use basechange_cublas !!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
            subroutine cumatmulnanoc_ZZ(A,devPtrX,C,M)
!---------------------------------------------------------------!
       integer , intent(in)  :: devPtrX
       integer, intent(in) :: M
       COMPLEX*16, intent(in) :: A(M,M)
       COMPLEX*16, intent(out) :: C(M,M)
       COMPLEX*16, dimension (:,:), ALLOCATABLE :: scratch
!---------------------------------------------------------------!
       allocate (scratch(M,M))
       call cumxp(A,devPtrX,scratch,M)
       call cumpxt(scratch,devPtrX,C,M)
       deallocate (scratch)
       return
       end
!===============================================================!
            subroutine cumatmulnanoc_CC(A,devPtrX,C,M)
!---------------------------------------------------------------!
       integer , intent(in)  :: devPtrX
       integer, intent(in) :: M
       COMPLEX*8, intent(in) :: A(M,M)
       COMPLEX*8, intent(out) :: C(M,M)
       COMPLEX*8, dimension (:,:), ALLOCATABLE :: scratch
!---------------------------------------------------------------!
       allocate (scratch(M,M))
       call cumxp(A,devPtrX,scratch,M)
       call cumpxt(scratch,devPtrX,C,M)
       deallocate (scratch)
       return
       end
!===============================================================!













