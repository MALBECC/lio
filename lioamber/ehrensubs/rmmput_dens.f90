!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine rmmput_dens_r(DensMao)
  use garcha_mod, only: Pmat_vec
  use basis_data, only: M
  implicit none
  real(kind=4),intent(in)     :: DensMao(M,M)
  integer               :: ii,jj,idx

  do jj=1,M
    do ii=jj,M
      idx=ii+(2*M-jj)*(jj-1)/2
      Pmat_vec(idx)=DensMao(ii,jj)*2
    enddo
    idx=jj+(2*M-jj)*(jj-1)/2
    Pmat_vec(idx)=Pmat_vec(idx)/2
  enddo

  return; end subroutine
!--------------------------------------------------------------------!
  subroutine rmmput_dens_d(DensMao)
  use garcha_mod, only: Pmat_vec
  use basis_data, only: M
  implicit none
  LIODBLE,intent(in)     :: DensMao(M,M)
  integer               :: ii,jj,idx

  do jj=1,M
    do ii=jj,M
      idx=ii+(2*M-jj)*(jj-1)/2
      Pmat_vec(idx)=DensMao(ii,jj)*2
    enddo
    idx=jj+(2*M-jj)*(jj-1)/2
    Pmat_vec(idx)=Pmat_vec(idx)/2
  enddo

  return; end subroutine
!--------------------------------------------------------------------!
  subroutine rmmput_dens_c(DensMao)
  use garcha_mod, only: Pmat_vec
  use basis_data, only: M
  implicit none
  complex(kind=4),intent(in)  :: DensMao(M,M)
  integer               :: ii,jj,idx

  do jj=1,M
    do ii=jj,M
      idx=ii+(2*M-jj)*(jj-1)/2
      Pmat_vec(idx)=(REAL(DensMao(ii,jj)))*2
    enddo
    idx=jj+(2*M-jj)*(jj-1)/2
    Pmat_vec(idx)=Pmat_vec(idx)/2
  enddo

  return; end subroutine
!--------------------------------------------------------------------!
  subroutine rmmput_dens_z(DensMao)
  use garcha_mod, only: Pmat_vec
  use basis_data, only: M
  implicit none
  complex(kind=8),intent(in) :: DensMao(M,M)
  integer               :: ii,jj,idx

  do jj=1,M
    do ii=jj,M
      idx=ii+(2*M-jj)*(jj-1)/2
      Pmat_vec(idx)=(REAL(DensMao(ii,jj)))*2
    enddo
    idx=jj+(2*M-jj)*(jj-1)/2
    Pmat_vec(idx)=Pmat_vec(idx)/2
  enddo

  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
