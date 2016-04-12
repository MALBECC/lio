!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine rmmget_fock_r(FockMao)
  use garcha_mod, only:M,RMM
  implicit none
  real*4,intent(out) :: FockMao(M,M)
  integer            :: ii,jj,idx,idx0

  idx0=M*(M+1)
  do jj=1,M
  do ii=jj,M
     idx=ii+(2*M-jj)*(jj-1)/2+idx0
     FockMao(ii,jj)=RMM(idx)
     FockMao(jj,ii)=RMM(idx)
  enddo
  enddo

  return; end subroutine
!--------------------------------------------------------------------!
  subroutine rmmget_fock_d(FockMao)
  use garcha_mod, only:M,RMM
  implicit none
  real*8,intent(out) :: FockMao(M,M)
  integer            :: ii,jj,idx,idx0

  idx0=M*(M+1)
  do jj=1,M
  do ii=jj,M
     idx=ii+(2*M-jj)*(jj-1)/2+idx0
     FockMao(ii,jj)=RMM(idx)
     FockMao(jj,ii)=RMM(idx)
  enddo
  enddo

  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
