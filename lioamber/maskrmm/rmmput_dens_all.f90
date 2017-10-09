!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine rmmput_dens_r(DensMao)
  use garcha_mod, only:M,RMM
  implicit none
  real*4,intent(in)     :: DensMao(M,M)
  integer               :: ii,jj,idx

  do jj=1,M
    do ii=jj,M
      idx=ii+(2*M-jj)*(jj-1)/2
      RMM(idx)=DensMao(ii,jj)*2
    enddo
    idx=jj+(2*M-jj)*(jj-1)/2
    RMM(idx)=RMM(idx)/2
  enddo

  return; end subroutine
!--------------------------------------------------------------------!
  subroutine rmmput_dens_d(DensMao)
  use garcha_mod, only:M,RMM
  implicit none
  real*8,intent(in)     :: DensMao(M,M)
  integer               :: ii,jj,idx

  do jj=1,M
    do ii=jj,M
      idx=ii+(2*M-jj)*(jj-1)/2
      RMM(idx)=DensMao(ii,jj)*2
    enddo
    idx=jj+(2*M-jj)*(jj-1)/2
    RMM(idx)=RMM(idx)/2
  enddo

  return; end subroutine
!--------------------------------------------------------------------!
  subroutine rmmput_dens_c(DensMao)
  use garcha_mod, only:M,RMM
  implicit none
  complex*8,intent(in)  :: DensMao(M,M)
  integer               :: ii,jj,idx

  do jj=1,M
    do ii=jj,M
      idx=ii+(2*M-jj)*(jj-1)/2
      RMM(idx)=(REAL(DensMao(ii,jj)))*2
    enddo
    idx=jj+(2*M-jj)*(jj-1)/2
    RMM(idx)=RMM(idx)/2
  enddo

  return; end subroutine
!--------------------------------------------------------------------!
  subroutine rmmput_dens_z(DensMao)
  use garcha_mod, only:M,RMM
  implicit none
  complex*16,intent(in) :: DensMao(M,M)
  integer               :: ii,jj,idx

  do jj=1,M
    do ii=jj,M
      idx=ii+(2*M-jj)*(jj-1)/2
      RMM(idx)=(REAL(DensMao(ii,jj)))*2
    enddo
    idx=jj+(2*M-jj)*(jj-1)/2
    RMM(idx)=RMM(idx)/2
  enddo

  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
