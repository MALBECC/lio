!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine rmmput_fock_r(FockMao)
  use garcha_mod, only: Fmat_vec
  use basis_data, only: M

  implicit none
  real*4,intent(in) :: FockMao(M,M)
  integer           :: ii,jj,idx

  do jj=1,M
  do ii=jj,M
     idx=ii+(2*M-jj)*(jj-1)/2
     Fmat_vec(idx) = FockMao(jj,ii)
  enddo
  enddo

  return; end subroutine
!--------------------------------------------------------------------!
  subroutine rmmput_fock_d(FockMao)
  use garcha_mod, only: Fmat_vec
  use basis_data, only: M

  implicit none
  real*8,intent(in) :: FockMao(M,M)
  integer           :: ii,jj,idx,idx0

  idx0=M*(M+1)
  do jj=1,M
  do ii=jj,M
     idx=ii+(2*M-jj)*(jj-1)/2+idx0
     Fmat_vec(idx) = FockMao(jj,ii)
  enddo
  enddo

  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
