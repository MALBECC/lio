!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine sdcmp_cholesky(Smat,Dvec,Vmat,Lmat,Linv,Umat,Uinv)
!--------------------------------------------------------------------!
       implicit none
       real*8,dimension(:,:),intent(in)  :: Smat
       real*8,dimension(:),intent(out)   :: Dvec
       real*8,dimension(:,:),intent(out) :: Vmat
       real*8,dimension(:,:),intent(out) :: Lmat,Linv
       real*8,dimension(:,:),intent(out) :: Umat,Uinv

       logical :: errorfound
       integer :: N,ii,jj,eid


! CHECK SIZES AND SYMETRY
!--------------------------------------------------------------------!
       N=size(Smat,1)
       errorfound=.false.
       errorfound=(errorfound).or.(size(Smat,2).ne.N)
       errorfound=(errorfound).or.(size(Dvec).ne.N)
       errorfound=(errorfound).or.(size(Vmat,1).ne.N)
       errorfound=(errorfound).or.(size(Vmat,2).ne.N)
       errorfound=(errorfound).or.(size(Lmat,1).ne.N)
       errorfound=(errorfound).or.(size(Lmat,2).ne.N)
       errorfound=(errorfound).or.(size(Linv,1).ne.N)
       errorfound=(errorfound).or.(size(Linv,2).ne.N)
       errorfound=(errorfound).or.(size(Umat,1).ne.N)
       errorfound=(errorfound).or.(size(Umat,2).ne.N)
       errorfound=(errorfound).or.(size(Uinv,1).ne.N)
       errorfound=(errorfound).or.(size(Uinv,2).ne.N)
       if (errorfound) stop('wrong sizes in cholesky_sdcmp')

       errorfound=.false.
       do jj=1,N
       do ii=1,jj
         errorfound=(errorfound).or.(Smat(ii,jj).ne.Smat(jj,ii))
       enddo
       enddo
       if (errorfound) stop('input matrix is not symetric')


! OBTAIN LMAT AND LINV - USING LAPACK
!--------------------------------------------------------------------!
!#ifdef pack
       Lmat=Smat
       call dpotrf('L', N, Lmat, N, eid)
       do ii=1,N-1
       do jj=ii+1,N
         Lmat(ii,jj)=0.0d0
       enddo
       enddo

       Linv=Lmat
       call dtrtri('L', 'N', N, Linv, N, eid)
!#endif


! OBTAIN EIGENVALUES AND EIGENVECTORS
!--------------------------------------------------------------------!
       do jj=1,N
         Dvec(jj)=Lmat(jj,jj)**2
         do ii=1,N
           Umat(jj,ii)=Lmat(ii,jj)
           Uinv(jj,ii)=Linv(ii,jj)
           Vmat(ii,jj)=Lmat(ii,jj)/Lmat(jj,jj)
         enddo
       enddo

       return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
