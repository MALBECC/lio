!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine rmmget_fock_r(FockMao)
  use garcha_mod, only: Fmat_vec
  use basis_data, only: M

  implicit none
  real(kind=4),intent(out) :: FockMao(M,M)
  integer            :: ii,jj,idx

  do jj=1,M
  do ii=jj,M
     idx=ii+(2*M-jj)*(jj-1)/2
     FockMao(ii,jj)=real(Fmat_vec(idx),4)
     FockMao(jj,ii)=real(Fmat_vec(idx),4)
  enddo
  enddo

  return; end subroutine
!--------------------------------------------------------------------!
  subroutine rmmget_fock_d(FockMao)
  use garcha_mod, only: Fmat_vec
  use basis_data, only: M

  implicit none
  LIODBLE,intent(out) :: FockMao(M,M)
  integer            :: ii,jj,idx

  do jj=1,M
  do ii=jj,M
     idx=ii+(2*M-jj)*(jj-1)/2
     FockMao(ii,jj) = Fmat_vec(idx)
     FockMao(jj,ii) = Fmat_vec(idx)
  enddo
  enddo

  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
