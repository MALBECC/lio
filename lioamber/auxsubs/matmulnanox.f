!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
! Hace C=Bt*(A*B) para matrices cuadradas
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine matmulnanox_d(Mati,Cmat,Mato,M)
       integer,intent(in)     :: M
       real*8,intent(in)      :: Cmat(M,M)
       real*8,intent(in)      :: Mati(M,M)
       real*8,intent(out)     :: Mato(M,M)
       real*8,allocatable     :: Matm(:,:),Trans(:,:)
       integer                :: ii,jj,kk

       allocate(Matm(M,M),Trans(M,M))
       do ii=1,M
       do jj=1,M
          Trans(ii,jj)=Mati(jj,ii)
       enddo
       enddo

       do ii=1,M
       do jj=1,M
         Matm(ii,jj)=DBLE(0)
         do kk=1,M
           Matm(ii,jj)=Matm(ii,jj)+Trans(kk,ii)*Cmat(kk,jj)
         enddo
       enddo
       enddo

       do ii=1,M
       do jj=1,M
         Mato(ii,jj)=DBLE(0)
         do kk=1,M
           Mato(ii,jj)=Mato(ii,jj)+Cmat(kk,ii)*Matm(kk,jj)
         enddo
       enddo
       enddo

       deallocate(Matm,Trans)

       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine matmulnanox_c(Mati,Cmat,Mato,M)
       integer,intent(in)     :: M
       real*8,intent(in)      :: Cmat(M,M)
       complex*8,intent(in)   :: Mati(M,M)
       complex*8,intent(out)  :: Mato(M,M)
       complex*8,allocatable  :: Matm(:,:),Trans(:,:)
       integer                :: ii,jj,kk

       allocate(Matm(M,M),Trans(M,M))
       do ii=1,M
       do jj=1,M
          Trans(ii,jj)=Mati(jj,ii)
       enddo
       enddo

       do ii=1,M
       do jj=1,M
         Matm(ii,jj)=CMPLX(0,0)
         do kk=1,M
           Matm(ii,jj)=Matm(ii,jj)+Trans(kk,ii)*Cmat(kk,jj)
         enddo
       enddo
       enddo

       do ii=1,M
       do jj=1,M
         Mato(ii,jj)=CMPLX(0,0)
         do kk=1,M
           Mato(ii,jj)=Mato(ii,jj)+Cmat(kk,ii)*Matm(kk,jj)
         enddo
       enddo
       enddo

       deallocate(Matm,Trans)

      return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
