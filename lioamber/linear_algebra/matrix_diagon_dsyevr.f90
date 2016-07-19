!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine matrix_diagon_dsyevr( matrix_in, eigen_vecs, eigen_vals , info )
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  real*8,  intent(in)            :: matrix_in(:,:)
  real*8,  intent(out)           :: eigen_vecs(:,:)
  real*8,  intent(out)           :: eigen_vals(:)
  integer, intent(out), optional :: info


  real*8, allocatable :: matrix_copy(:,:)
  integer,allocatable :: neig_vecs(:)
  integer             :: neig_vals
  integer             :: lwork
  real*8, allocatable :: work(:)
  integer             :: liwork
  integer,allocatable :: iwork(:)

  integer            :: M
  integer            :: local_stat
  real*8             :: dlamch
  real*8,  parameter :: zero_d=0.0d0
  integer, parameter :: zero_i=0
!
!
!
!
! Initial checks
!------------------------------------------------------------------------------!
  local_stat=0
  if ( present(info) ) info=0

  M=size(matrix_in,1)
  if ( M /= size(matrix_in,2) )  local_stat=1
  if ( M /= size(eigen_vecs,1) ) local_stat=2
  if ( M /= size(eigen_vecs,2) ) local_stat=3
  if ( M /= size(eigen_vals) )   local_stat=4

  if ( local_stat /= 0 ) then
    if ( present(info) ) then
      info = 1
      return
    else
      print*,'matrix_diagon_dsyevr : incompatible size between arguments'
      print*,'local info: ', local_stat
      stop
    end if
  end if
!
!
!
!
! Allocation of variables
!------------------------------------------------------------------------------!
  if ( local_stat == 0 ) then
    if ( allocated(work) ) deallocate(work)
    allocate( work(1), stat=local_stat )
  end if

  if ( local_stat == 0 ) then
    if ( allocated(iwork) ) deallocate(iwork)
    allocate( iwork(1), stat=local_stat )
  end if

  if ( local_stat == 0 ) then
    if ( allocated(neig_vecs) ) deallocate(neig_vecs)
    allocate( neig_vecs(2*M), stat=local_stat )
  end if

  if ( local_stat == 0 ) then
    if ( allocated(matrix_copy) ) deallocate(matrix_copy)
    allocate( matrix_copy(M,M), stat=local_stat )
  end if

  if ( local_stat /= 0 ) then
    if ( present(info) ) then
      info = 2
      return
    else
      print*,'matrix_diagon_dsyevr : critical error during initial allocation'
      print*,'local info: ', local_stat
      stop
    end if
  end if
!
!
!
!
! Set working size
!------------------------------------------------------------------------------!
  lwork=-1
  liwork=-1
  matrix_copy=matrix_in

  call dsyevr( 'V', 'A', 'L', M, matrix_copy, M, zero_d, zero_d, &
               zero_i, zero_i, dlamch('Safe minimum'), neig_vals, &
               eigen_vals, eigen_vecs, M, neig_vecs, work, lwork, &
               iwork, liwork, local_stat)


  if ( local_stat == 0 ) then
    lwork = work(1)
    if ( allocated(work) ) deallocate( work )
    allocate( work(lwork), stat=local_stat )
  end if

  if ( local_stat == 0 ) then
    liwork = iwork(1)
     if ( allocated(iwork) ) deallocate( iwork )
    allocate( iwork(liwork), stat=local_stat  )
  end if

  if ( local_stat /= 0 ) then
    if ( present(info) ) then
      info = 3
      return
    else
      print*,'matrix_diagon_dsyevr : critical error while setting working size'
      print*,'local info: ', local_stat
      stop
    end if
  end if
!
!
!
!
! Do actual diagonalization
!------------------------------------------------------------------------------!
  call dsyevr( 'V', 'A', 'L', M, matrix_copy, M, zero_d, zero_d, &
               zero_i, zero_i, dlamch('Safe minimum'), neig_vals, &
               eigen_vals, eigen_vecs, M, neig_vecs, work, lwork, &
               iwork, liwork, local_stat)

  if ( local_stat /= 0 ) then
    if ( present(info) ) then
      info = 4
      return
    else
      print*,'matrix_diagon_dsyevr : critical error while diagonalizing'
      print*,'local info: ', local_stat
      stop
    end if
  end if
!
!
!
!
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
