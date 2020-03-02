!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine matrix_diagon_dsyevd( matrix_in, eigen_vecs, eigen_vals , info )
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  LIODBLE,  intent(in)            :: matrix_in(:,:)
  LIODBLE,  intent(out)           :: eigen_vecs(:,:)
  LIODBLE,  intent(out)           :: eigen_vals(:)
  integer, intent(out), optional :: info

  integer             :: lwork
  LIODBLE, allocatable :: work(:)
  integer             :: liwork
  integer,allocatable :: iwork(:)

  integer :: M,ii,jj
  integer :: local_stat
  logical :: be_safe
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

  if ( local_stat /= 0 ) then
    if ( present(info) ) then
      info = 1
      return
    else
      print*,'matrix_diagon_dsyevd : critical error during initial allocation'
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
  eigen_vecs=matrix_in

# ifdef magma
  call magmaf_dsyevd( 'V', 'L', M, eigen_vecs, M, eigen_vals, &
                    & work, lwork, iwork, lwork, local_stat )
# else
  call        dsyevd( 'V', 'L', M, eigen_vecs, M, eigen_vals, &
                    & work, lwork, iwork, lwork, local_stat )
# endif

  if ( local_stat == 0 ) then
    lwork = int(work(1))
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
      info = 2
      return
    else
      print*,'matrix_diagon_dsyevd : critical error while setting working size'
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
# ifdef magma
  call magmaf_dsyevd( 'V', 'L', M, eigen_vecs, M, eigen_vals, &
                    & work, lwork, iwork, liwork, local_stat )
# else
  call        dsyevd( 'V', 'L', M, eigen_vecs, M, eigen_vals, &
                    & work, lwork, iwork, liwork, local_stat )
# endif

  if ( local_stat /= 0 ) then
    if ( present(info) ) then
      info = 3
      return
    else
      print*,'matrix_diagon_dsyevd : critical error while diagonalizing'
      print*,'local info: ', local_stat
      stop
    end if
  end if

  be_safe=.true.
  if (be_safe.eqv..true.) then
    do ii=1,size(eigen_vecs,1)
    do jj=1,size(eigen_vecs,2)
      if ( eigen_vecs(ii,jj) /= eigen_vecs(ii,jj) ) then
        if ( present(info) ) then
          info = 4
          return
        else
          print*,'matrix_diagon_dsyevd : critical error has not been &
                 &recognized by lapack!'
          print*,'local info: ', local_stat
          stop
        end if
      end if
    end do
    end do
  end if
!
!
!
!
end subroutine matrix_diagon_dsyevd
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
