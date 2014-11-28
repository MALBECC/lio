!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       program testprog
       use auxsubs
       implicit none
       real*8,allocatable     :: realvec(:)
       complex*16,allocatable :: compvec(:)
       real*8,allocatable     :: realmat(:,:)
       complex*16,allocatable :: compmat(:,:)
       real*8                 :: pos,mid,a

       allocate(realvec(3),compvec(3))
       realvec(1)=DBLE(1.0)
       realvec(2)=DBLE(1.0)
       realvec(3)=DBLE(1.0)
       compvec(1)=DCMPLX(1.0,2.0)
       compvec(2)=DCMPLX(1.0,2.0)
       compvec(3)=DCMPLX(1.0,2.0)

       pos=1.0d0
       mid=0.0d0
       a=1.0d0
       print*,gaussbell(pos,mid,a,realvec)
       print*,gaussbell(pos,mid,a,compvec)

       allocate(realmat(2,2),compmat(2,2))
       realmat(1,1)=DBLE(0.0)
       realmat(1,2)=DBLE(1.0)
       realmat(2,1)=DBLE(1.0)
       realmat(2,2)=DBLE(0.0)
       compmat(1,1)=DCMPLX(2.0,0.0)
       compmat(1,2)=DCMPLX(0.0,0.0)
       compmat(2,1)=DCMPLX(0.0,0.0)
       compmat(2,2)=DCMPLX(1.0,0.0)

       print*,commutate(realmat,compmat)
       print*,commutate(realmat,realmat)


       stop;end program
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
