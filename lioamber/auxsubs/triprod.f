!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function triprod_ddd(MIzq,MCen,MDer)
     > result(MOut)
       implicit none
       real*8,intent(in)      :: MIzq(:,:),MDer(:,:)
       real*8,intent(in)      :: MCen(:,:)
       real*8,allocatable     :: MOut(:,:)
       real*8                 :: pc
       real*8                 :: pi,pd
       integer                :: nii,nni,nnd,ndd
     > ,                         kii,kki,kkd,kdd

       nii=size(MIzq,1)
       nni=size(MIzq,2)
       nnd=size(MDer,1)
       ndd=size(MDer,2)
       allocate(MOut(nii,ndd))

       do kii=1,nii
       do kdd=1,ndd
         MOut(kii,kdd)=0.0d0
         do kki=1,nni
         do kkd=1,nnd
           pi=MIzq(kii,kki)
           pc=MCen(kki,kkd)
           pd=MDer(kkd,kdd)
           MOut(kii,kdd)=MOut(kii,kdd)+pi*pc*pd
         enddo
         enddo
       enddo
       enddo
       return;end function
!
!
!--------------------------------------------------------------------!
       function triprod_dcd(MIzq,MCen,MDer)
     > result(MOut)
       implicit none
       real*8,intent(in)      :: MIzq(:,:),MDer(:,:)
       complex*16,intent(in)  :: MCen(:,:)
       complex*16,allocatable :: MOut(:,:)
       complex*16             :: pc
       real*8                 :: pi,pd
       integer                :: nii,nni,nnd,ndd
     > ,                         kii,kki,kkd,kdd

       nii=size(MIzq,1)
       nni=size(MIzq,2)
       nnd=size(MDer,1)
       ndd=size(MDer,2)
       allocate(MOut(nii,ndd))

       do kii=1,nii
       do kdd=1,ndd
         MOut(kii,kdd)=DCMPLX(0.0d0,0.0d0)
         do kki=1,nni
         do kkd=1,nnd
           pi=MIzq(kii,kki)
           pc=MCen(kki,kkd)
           pd=MDer(kkd,kdd)
           MOut(kii,kdd)=MOut(kii,kdd)+pi*pc*pd
         enddo
         enddo
       enddo
       enddo
       return;end function
!
!
!--------------------------------------------------------------------!
       function triprod_ccc(MIzq,MCen,MDer)
     > result(MOut)
       implicit none
       complex*16,intent(in)  :: MIzq(:,:),MDer(:,:)
       complex*16,intent(in)  :: MCen(:,:)
       complex*16,allocatable :: MOut(:,:)
       complex*16             :: pc
       complex*16             :: pi,pd
       integer                :: nii,nni,nnd,ndd
     > ,                         kii,kki,kkd,kdd

       nii=size(MIzq,1)
       nni=size(MIzq,2)
       nnd=size(MDer,1)
       ndd=size(MDer,2)
       allocate(MOut(nii,ndd))

       do kii=1,nii
       do kdd=1,ndd
         MOut(kii,kdd)=DCMPLX(0.0d0,0.0d0)
         do kki=1,nni
         do kkd=1,nnd
           pi=MIzq(kii,kki)
           pc=MCen(kki,kkd)
           pd=MDer(kkd,kdd)
           MOut(kii,kdd)=MOut(kii,kdd)+pi*pc*pd
         enddo
         enddo
       enddo
       enddo
       return;end function
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
