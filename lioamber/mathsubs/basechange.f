!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! BASECHANGE
! BASETRANSFORM PROCEDURES
!
! (1) Initialization of Matm(nnd,ndd) and Mato(nii,ndd)
! (2) First Product Mati(nni,nnd)*Utrp(nnd,ndd)
! (3) Second Product Utrp(nii,nni)*Matm(nni,ndd)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function basechange_d(M,Utrp,Mati) result(Mato)
       implicit none
       integer,intent(in)     :: M
       real*8,intent(in)      :: Utrp(M,M)
       real*8,intent(in)      :: Mati(M,M)
       real*8,allocatable     :: Matm(:,:)
       real*8,allocatable     :: Mato(:,:)
       integer                :: kii,kki,kkd,kdd
!
       allocate(Matm(M,M),Mato(M,M))
       Matm=DBLE(0)
       Mato=DBLE(0)
!
       do kdd=1,M
       do kkd=1,M
       do kki=1,M
         Matm(kki,kdd)=Matm(kki,kdd)+Mati(kki,kkd)*Utrp(kdd,kkd)
       enddo
       enddo
       enddo
!
       do kdd=1,M
       do kki=1,M
       do kii=1,M
         Mato(kii,kdd)=Mato(kii,kdd)+Utrp(kii,kki)*Matm(kki,kdd)
       enddo
       enddo
       enddo
!
       return;end function
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function basechange_c(M,Utrp,Mati) result(Mato)
       implicit none
       integer,intent(in)     :: M
       real*8,intent(in)      :: Utrp(M,M)
       complex*8,intent(in)   :: Mati(M,M)
       complex*8,allocatable  :: Matm(:,:)
       complex*8,allocatable  :: Mato(:,:)
       integer                :: kii,kki,kkd,kdd
!
       allocate(Matm(M,M),Mato(M,M))
       Matm=DCMPLX(0,0)
       Mato=DCMPLX(0,0)
!
       do kdd=1,M
       do kkd=1,M
       do kki=1,M
         Matm(kki,kdd)=Matm(kki,kdd)+Mati(kki,kkd)*Utrp(kdd,kkd)
       enddo
       enddo
       enddo
!
       do kdd=1,M
       do kki=1,M
       do kii=1,M
         Mato(kii,kdd)=Mato(kii,kdd)+Utrp(kii,kki)*Matm(kki,kdd)
       enddo
       enddo
       enddo
!
       return;end function
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function basechange_z(M,Utrp,Mati) result(Mato)
       implicit none
       integer,intent(in)     :: M
       real*8,intent(in)      :: Utrp(M,M)
       complex*16,intent(in)  :: Mati(M,M)
       complex*16,allocatable :: Matm(:,:)
       complex*16,allocatable :: Mato(:,:)
       integer                :: kii,kki,kkd,kdd
!
       allocate(Matm(M,M),Mato(M,M))
       Matm=DCMPLX(0,0)
       Mato=DCMPLX(0,0)
!
       do kdd=1,M
       do kkd=1,M
       do kki=1,M
         Matm(kki,kdd)=Matm(kki,kdd)+Mati(kki,kkd)*Utrp(kdd,kkd)
       enddo
       enddo
       enddo
!
       do kdd=1,M
       do kki=1,M
       do kii=1,M
         Mato(kii,kdd)=Mato(kii,kdd)+Utrp(kii,kki)*Matm(kki,kdd)
       enddo
       enddo
       enddo
!
       return;end function
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
