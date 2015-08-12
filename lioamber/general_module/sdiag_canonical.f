!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine sdiag_canonical(Smat,Dvec,Vmat,Xmat,Xtrp,Ymat,Ytrp)
!--------------------------------------------------------------------!
       implicit none
       real*8,dimension(:,:),intent(in)  :: Smat
       real*8,dimension(:),intent(out)   :: Dvec
       real*8,dimension(:,:),intent(out) :: Vmat
       real*8,dimension(:,:),intent(out) :: Xmat,Xtrp
       real*8,dimension(:,:),intent(out) :: Ymat,Ytrp

       logical :: errorfound
       integer :: N,ii,jj,eid
       real*8,parameter :: zero_value=1.0d-06

       real*8,dimension(:),allocatable :: WORK
       integer                         :: LWORK



! CHECK SIZES AND SYMETRY
!--------------------------------------------------------------------!
       N=size(Smat,1)
       errorfound=.false.
       errorfound=(errorfound).or.(size(Smat,2).ne.N)
       errorfound=(errorfound).or.(size(Dvec).ne.N)
       errorfound=(errorfound).or.(size(Vmat,1).ne.N)
       errorfound=(errorfound).or.(size(Vmat,2).ne.N)
       errorfound=(errorfound).or.(size(Xmat,1).ne.N)
       errorfound=(errorfound).or.(size(Xmat,2).ne.N)
       errorfound=(errorfound).or.(size(Xtrp,1).ne.N)
       errorfound=(errorfound).or.(size(Xtrp,2).ne.N)
       errorfound=(errorfound).or.(size(Ymat,1).ne.N)
       errorfound=(errorfound).or.(size(Ymat,2).ne.N)
       errorfound=(errorfound).or.(size(Ytrp,1).ne.N)
       errorfound=(errorfound).or.(size(Ytrp,2).ne.N)
       if (errorfound) stop('wrong sizes in canonical_sdiag')

       errorfound=.false.
       do jj=1,N
       do ii=1,jj
         errorfound=(errorfound).or.(Smat(ii,jj).ne.Smat(jj,ii))
       enddo
       enddo
       if (errorfound) stop('input matrix is not symetric')


! DIAGONALIZE SMAT - USING LAPACK
!--------------------------------------------------------------------!
!#ifdef pack
       Vmat=Smat
       allocate(WORK(1))
       call dsyev('V', 'L', N, Vmat, N, Dvec, WORK, -1, eid)
       LWORK=int(WORK(1))
       deallocate(WORK)
       allocate(WORK(LWORK))
       call dsyev('V', 'L', N, Vmat, N, Dvec, WORK, LWORK, eid)
!#endif


! CONSTRUCT XMAT AND YMAT
!--------------------------------------------------------------------!
       do jj=1,N
         if (Dvec(jj).lt.zero_value) then
           write(*,*) 'LINEAR DEPENDENCY DETECTED'
           do ii=1,N
             Xmat(ii,jj)=0.0d0
             Xtrp(jj,ii)=0.0d0
             Ymat(ii,jj)=0.0d0
             Ytrp(jj,ii)=0.0d0
           enddo
         else
           do ii=1,N
             Xmat(ii,jj)=Vmat(ii,jj)/sqrt(Dvec(jj))
             Xtrp(jj,ii)=Xmat(ii,jj)
             Ymat(ii,jj)=Xmat(ii,jj)*Dvec(jj)
             Ytrp(jj,ii)=Ymat(ii,jj)
           enddo
         endif
       enddo

       return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
