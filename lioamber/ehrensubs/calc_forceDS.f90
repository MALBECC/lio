!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine calc_forceDS &
  ( Natoms, Nbasis, nucpos, nucvel, DensMao, FockMao, Sinv, Bmat, forceDS )
!------------------------------------------------------------------------------!
!
! DESCRIPTION
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in)    :: Natoms
  integer,intent(in)    :: Nbasis
  real*8,intent(in)     :: nucpos(3,Natoms)
  real*8,intent(in)     :: nucvel(3,Natoms)
  complex*16,intent(in) :: DensMao(Nbasis,Nbasis)
  real*8,intent(in)     :: FockMao(Nbasis,Nbasis)
  real*8,intent(in)     :: Sinv(Nbasis,Nbasis)
  real*8,intent(out)    :: Bmat(Nbasis,Nbasis)
  real*8,intent(out)    :: forceDS(3,Natoms)


  real*8,allocatable     :: Btrp(:,:)
  complex*16,allocatable :: InputMat(:,:),MatTrp(:,:),MatDir(:,:)
  complex*16,allocatable :: fterm1(:,:),fterm2(:,:),fterm3(:,:)
!
!
!------------------------------------------------------------------------------!
  call g2g_timer_start('calc_forceDS')
  allocate(InputMat(Nbasis,Nbasis))
  allocate(MatTrp(Nbasis,Nbasis),MatDir(Nbasis,Nbasis))
  allocate(Btrp(Nbasis,Nbasis))
  allocate(fterm1(3,Natoms),fterm2(3,Natoms),fterm3(3,Natoms))


  fterm1=dcmplx(0.0d0,0.0d0)
  fterm2=dcmplx(0.0d0,0.0d0)
  fterm3=dcmplx(0.0d0,0.0d0)

! NOTA: El orden de las multiplicaciones afecta levemente el
! resultado obtenido
  MatTrp=matmul(DensMao,FockMao)
  MatTrp=matmul(MatTrp,Sinv)
  MatDir=matmul(FockMao,DensMao)
  MatDir=matmul(Sinv,MatDir)
  InputMat=transpose(MatTrp)+MatDir
  call calc_forceDS_dss(Natoms,Nbasis,nucpos,nucvel,InputMat,Bmat,fterm1)
  Btrp=transpose(Bmat)


  MatTrp=matmul(DensMao,Btrp)
  MatTrp=matmul(MatTrp,Sinv)
  MatTrp=MatTrp*dcmplx(0.0d0, 1.0d0)
  MatDir=matmul(Sinv,Bmat)
  MatDir=matmul(MatDir,DensMao)
  MatDir=MatDir*dcmplx(0.0d0,-1.0d0)
  InputMat=transpose(MatTrp)+MatDir
  call calc_forceDS_dss(Natoms,Nbasis,nucpos,nucvel,InputMat,Bmat,fterm2)


  MatTrp=DensMao*dcmplx(0.0d0,-1.0d0)
  MatDir=DensMao*dcmplx(0.0d0, 1.0d0)
  InputMat=transpose(MatTrp)+MatDir
  call calc_forceDS_dds(Natoms,Nbasis,nucpos,nucvel,InputMat,fterm3)


  forceDS=dble(real(fterm1+fterm2+fterm3))

  deallocate(InputMat,MatTrp,MatDir,Btrp)
  deallocate(fterm1,fterm2,fterm3)
  call g2g_timer_stop('calc_forceDS')
end subroutine calc_forceDS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
