!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine testforce(Sinv,FockMao,DensMao)
!--------------------------------------------------------------------!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use garcha_mod, &
      only:M,a,c,r,nshell,natom,ncont,nl,nuc,RMM,Md,nucvel
  use ehrenfest
  use basis_data, &
      only:basis_data_set
  use testmod

  implicit none
  real*8,intent(in)      :: Sinv(M,M)
  real*8,intent(in)      :: FockMao(M,M)
  complex*16,intent(in)  :: DensMao(M,M)

  real*8,allocatable     :: nr(:,:),nv(:,:),as(:,:),cs(:,:)
  real*8,allocatable     :: ffold(:,:),ffnew(:,:),Mato(:,:)
  complex*16,allocatable :: AuxMat(:,:), Mdir(:,:),Mtrp(:,:)
  complex*16,allocatable :: Bmat(:,:),ffalt(:,:)
  integer,allocatable    :: nucof(:)
  integer                :: kk,ii,jj,unitid,nk,nn
  integer                :: MM,MMd,indx

!--------------------------------------------------------------------!
  allocate(nr(3,natom),nv(3,natom),as(nl,M),cs(nl,M),nucof(M))
  allocate(AuxMat(M,M),Mdir(M,M),Mtrp(M,M),Bmat(M,M))
  allocate(ffold(natom,3),ffnew(3,natom),ffalt(3,natom))
  allocate(Mato(M,M))

  call testinit()
  call g2g_timer_start('oldie')
  DSX(:,:,:)=0.0d0
  DSY(:,:,:)=0.0d0
  DSZ(:,:,:)=0.0d0
  ffold=0.0d0
  call calcDSM(ffold)
  call g2g_timer_stop('oldie')
!  call intSG(ffold)
  do kk=1,natom
   do ii=1,M
   do jj=1,M
     DSXc(ii,jj,kk)=DSX(ii,jj,kk)
     DSYc(ii,jj,kk)=DSY(ii,jj,kk)
     DSZc(ii,jj,kk)=DSZ(ii,jj,kk)
   enddo
   enddo
  enddo


! Crear Mdir y Mtrp
  Mdir=matmul(FockMao,DensMao)
  Mdir=matmul(Sinv,Mdir)
  Mato=(-1)*DBLE(Mdir)
  Mtrp=matmul(DensMao,FockMao)
  Mtrp=matmul(Mtrp,Sinv)
  Mtrp=transpose(Mtrp)
  AuxMat=Mdir+Mtrp

! Transpone a y c
  do ii=1,M
    nucof(ii)=nuc(ii)
    do jj=1,nl
      as(jj,ii)=a(ii,jj)
      cs(jj,ii)=c(ii,jj)
    enddo
  enddo

! Transpone r y setea nv=0
  do ii=1,natom
  do jj=1,3
    nr(jj,ii)=r(ii,jj)
    nv(jj,ii)=0.0d0
  enddo
  enddo

! Calcula las fuerzas
  DSX(:,:,:)=0.0d0
  DSY(:,:,:)=0.0d0
  DSZ(:,:,:)=0.0d0
  ffnew=DCMPLX(0.0d0,0.0d0)
  call g2g_timer_start('newbie')
  call fzaDS2(natom,M,nshell(0),nshell(1),ncont,nl, &
              AuxMat,nr,nv,as,cs,nucof,Bmat,ffalt)
  call g2g_timer_stop('newbie')
  do kk=1,natom
   do ii=1,M
   do jj=1,M
     DSXt(ii,jj,kk)=DSX(jj,ii,kk)
     DSYt(ii,jj,kk)=DSY(jj,ii,kk)
     DSZt(ii,jj,kk)=DSZ(jj,ii,kk)
   enddo
   enddo
  enddo

  DSX=DSX+DSXt
  DSY=DSY+DSYt
  DSZ=DSZ+DSZt


  print*,'     natom    dir'
  do ii=1,natom
  do kk=1,3
    print*,ii,kk,ffold(ii,kk),ffalt(kk,ii)
  enddo
  enddo

  print*,''
  print*,''
  call basis_data_set&
       (nshell(0),nshell(1),nshell(2),nucof,ncont,a,c)
  call calc_forceDS(natom,M,nr,nucvel,DensMao,FockMao,Sinv,Mato,ffnew)
  print*,'     natom    dir'
  do ii=1,natom
  do kk=1,3
    print*,ii,kk,ffold(ii,kk),ffnew(kk,ii)
  enddo
  enddo
  print*,''
  print*,''


!--------------------------------------------------------------------!
  deallocate(Mato)
  deallocate(nr,nv,as,cs,nucof)
  deallocate(AuxMat,Mdir,Mtrp,Bmat)
  deallocate(ffold,ffnew)
  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
