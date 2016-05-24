!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fockon_diagonalize(M,fock_on,eig_vecs,eig_vals)
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer, intent(in)  :: M
  real*8,  intent(in)  :: fock_on(M,M)
  real*8,  intent(out) :: eig_vecs(M,M)
  real*8,  intent(out) :: eig_vals(M)

  real*8, allocatable :: fock_aux(:,:)
  integer,allocatable :: neig_vecs(:)
  integer             :: neig_vals
  integer             :: lwork, liwork, info
  real*8, allocatable :: work(:)
  integer,allocatable :: iwork(:)

  real*8             :: dlamch
  real*8,  parameter :: zero_d=0.0d0
  integer, parameter :: zero_i=0
!
!
!
! Setup for work and iwork
!------------------------------------------------------------------------------!
  if ( allocated(fock_aux) )  deallocate( fock_aux )
  allocate( fock_aux(M,M) )
  fock_aux=fock_on

  if ( allocated(neig_vecs) )  deallocate( neig_vecs )
  allocate( neig_vecs(2*M) )

  if ( allocated(work) )  deallocate( work );  allocate( work(1) )
  if ( allocated(iwork) ) deallocate( iwork ); allocate( iwork(1) )
  info=0
  lwork=-1
  liwork=-1
  eig_vecs=fock_on

# ifdef magma
  call magmaf_dsyevd( 'V', 'L', M, eig_vecs, M, eig_vals, &
       work, lwork, iwork, lwork, info )
# else
!  call dsyevd( 'V', 'L', M, eig_vecs, M, eig_vals, &
!       work, lwork, iwork, lwork, info )
  call dsyevr( 'V', 'A', 'L', M, fock_aux, M, zero_d, zero_d, zero_i, zero_i, &
       dlamch('Safe minimum'), neig_vals, eig_vals, eig_vecs, M, neig_vecs, &
       work, lwork, iwork, liwork, info)
# endif

  if ( info .ne. 0 ) then
    print*,'Error while diagonalizing fock (pre): ',info
    stop
  endif

  lwork = work(1)
  deallocate( work )
  allocate( work(lwork) )

  liwork = iwork(1)
  deallocate( iwork )
  allocate( iwork(liwork) )
!
!
! Actual diagonalization
!------------------------------------------------------------------------------!
# ifdef magma
  call magmaf_dsyevd( 'V', 'L', M, eig_vecs, M, eig_vals, &
       work, lwork, iwork, liwork, info)
# else
!  call dsyevd( 'V', 'L', M, eig_vecs, M, eig_vals, &
!       work, lwork, iwork, liwork, info )
  call dsyevr( 'V', 'A', 'L', M, fock_aux, M, zero_d, zero_d, zero_i, zero_i, &
       dlamch('Safe minimum'), neig_vals, eig_vals, eig_vecs, M, neig_vecs, &
       work, lwork, iwork, liwork, info)
# endif

  if ( info /= 0 ) then
    print*,'Error while diagonalizing fock (post): ',info
    stop
  endif

end subroutine fockon_diagonalize
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
