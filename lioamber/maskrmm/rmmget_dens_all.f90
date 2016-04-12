!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine rmmget_dens_r(DensMao)
  use garcha_mod, only:M,RMM
  implicit none
  real*4,intent(out)     :: DensMao(M,M)
  integer                :: ii,jj,idx

  do jj=1,M
    do ii=jj,M
      idx=ii+(2*M-jj)*(jj-1)/2
      DensMao(ii,jj)=RMM(idx)
      DensMao(jj,ii)=RMM(idx)
    enddo
    DensMao(jj,jj)=DensMao(jj,jj)/2.0d0
  enddo

  return;end subroutine
!--------------------------------------------------------------------!
  subroutine rmmget_dens_d(DensMao)
  use garcha_mod, only:M,RMM
  implicit none
  real*8,intent(out)     :: DensMao(M,M)
  integer                :: ii,jj,idx

  do jj=1,M
    do ii=jj,M
      idx=ii+(2*M-jj)*(jj-1)/2
      DensMao(ii,jj)=RMM(idx)
      DensMao(jj,ii)=RMM(idx)
    enddo
    DensMao(jj,jj)=DensMao(jj,jj)/2.0d0
  enddo

  return;end subroutine
!--------------------------------------------------------------------!
  subroutine rmmget_dens_c(DensMao)
  use garcha_mod, only:M,RMM
  implicit none
  complex*8,intent(out)  :: DensMao(M,M)
  integer                :: ii,jj,idx

  do jj=1,M
    do ii=jj,M
      idx=ii+(2*M-jj)*(jj-1)/2
      DensMao(ii,jj)=cmplx(RMM(idx))
      DensMao(jj,ii)=cmplx(RMM(idx))
    enddo
    DensMao(jj,jj)=DensMao(jj,jj)/2.0d0
  enddo

  return;end subroutine
!--------------------------------------------------------------------!
  subroutine rmmget_dens_z(DensMao)
  use garcha_mod, only:M,RMM
  implicit none
  complex*16,intent(out) :: DensMao(M,M)
  integer                :: ii,jj,idx

  do jj=1,M
    do ii=jj,M
      idx=ii+(2*M-jj)*(jj-1)/2
      DensMao(ii,jj)=dcmplx(RMM(idx))
      DensMao(jj,ii)=dcmplx(RMM(idx))
    enddo
    DensMao(jj,jj)=DensMao(jj,jj)/2.0d0
  enddo

  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
