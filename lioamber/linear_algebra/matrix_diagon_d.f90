!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine matrix_diagon_d( matrix_in, eigen_vecs, eigen_vals , info )
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  real*8,  intent(in)            :: matrix_in(:,:)
  real*8,  intent(out)           :: eigen_vecs(:,:)
  real*8,  intent(out)           :: eigen_vals(:)
  integer, intent(out), optional :: info
  integer                        :: local_stat

  local_stat=0
  if ( present(info) ) info=0

! WARNING: I had cases where dsyevd fails but no information is passed to
!          local stat. May be necessary to always use dsyevr...
  call matrix_diagon_dsyevd &
       ( matrix_in, eigen_vecs, eigen_vals , info=local_stat )

  if ( local_stat /= 0 ) then
    if ( present(info) ) then
      info = -1
    end if
    local_stat=0
    call matrix_diagon_dsyevr &
         ( matrix_in, eigen_vecs, eigen_vals , info=local_stat )
  end if


  if ( local_stat /= 0 ) then
    if ( present(info) ) then
      info = 1
      return
    else
      print*,'matrix_diagon_d : all diagonalization attempts failed'
      print*,'local info: ', local_stat
      stop
    end if
  end if

end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
