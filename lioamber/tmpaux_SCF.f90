!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module tmpaux_SCF
   implicit none
   contains
!
! TODO:
!
!  neighbor_list_2e: Find out what this is about.
!
!  starting_guess: Generates density starting guess. Find out how it build
!                  and add a better description.
!
!  COPY_VEC: this should go inside of maskrmm, or replaced by a subroutine
!            there,
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   SUBROUTINE neighbor_list_2e()
! Para hacer lineal la integral de 2 electrone con lista de vecinos. Nano
      USE garcha_mod, ONLY : natom, natomc, r, d, jatc, rmax, nshell, atmin,   &
                           & nnps, nnpp, nnpd, M, nuc
      IMPLICIT NONE
      INTEGER :: i,j, iij, iik, iikk
      REAL*8  :: zij, ti, tj, alf, rexp
      
      do i=1,natom
         natomc(i)=0
         do j=1,natom
            d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+(r(i,3)-r(j,3))**2
            zij=atmin(i)+atmin(j)
            ti=atmin(i)/zij
            tj=atmin(j)/zij
            alf=atmin(i)*tj
            rexp=alf*d(i,j)
            if (rexp.lt.rmax) then
               natomc(i)=natomc(i)+1
               jatc(natomc(i),i)=j
            endif
         enddo
      enddo

      do iij=nshell(0),1,-1
        nnps(nuc(iij))=iij
      enddo

      do iik=nshell(0)+nshell(1),nshell(0)+1,-1
        nnpp(nuc(iik))=iik
      enddo

      do iikk=M,nshell(0)+nshell(1)+1,-1
        nnpd(nuc(iikk))=iikk
      enddo
   END SUBROUTINE neighbor_list_2e



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   SUBROUTINE starting_guess_old(xnano)
      use garcha_mod, ONLY: RMM, VCINP, primera, M, X, Md, NCO, MO_coef_at
      IMPLICIT NONE
      integer :: info
      real*8, dimension (M,M), intent(inout)::xnano
      real*8, dimension (:), ALLOCATABLE :: rmm5,rmm15
      integer :: M1,M2,M3, M5, M7, M9, M11, M13, M15, M17, M18, MM, MMd 
               ! temporales hasta q rompamos RMM
      integer :: ii,jj,i,j,k,kk !auxiliares
      real*8 :: ff

      call g2g_timer_start('initial guess')
      call g2g_timer_sum_stop('Overlap decomposition')

      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2
      allocate(rmm5(MM),rmm15(mm))

      M1=1 ! first P
      M2=2*M
      M3=M1+MM ! now Pnew
      M5=M3+MM! now S, F also uses the same position after S was used
      M7=M5+MM! now G
      M9=M7+MMd ! now Gm
      M11=M9+MMd! now H
      M13=M11+MM! W ( eigenvalues ), also this space is used in least squares
      M15=M13+M! aux ( vector for ESSl)
      M17=M15+MM! Least squares
      M18=M17+MMd! vectors of MO

!     CASE OF NO STARTING GUESS PROVIDED, 1 E FOCK MATRIX USED
!     FCe = SCe; (X^T)SX = 1
!     F' = (X^T)FX
!     => (X^-1*C)^-1 * F' * (X^-1*C) = e
!
!     Calculate F' in RMM(M5)
      if ((.not.VCINP).and.primera) then
         call g2g_timer_sum_start('initial guess')
         primera=.false.
         do i=1,M
         do j=1,M
            X(i,M+j)=0.D0
            do k=1,j
               X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M11+j+(M2-k)*(k-1)/2-1)
            enddo
            do k=j+1,M
               X(i,M+j)=X(i,M+j)+X(k,i)*RMM(M11+k+(M2-j)*(j-1)/2-1)
            enddo
         enddo
         enddo


         kk=0
         do j=1,M
         do i=j,M
            kk=kk+1
            RMM(M5+kk-1)=0.D0
            do k=1,j
              RMM(M5+kk-1)=RMM(M5+kk-1)+X(i,M+k)*X(k,j)
            enddo
         enddo
         enddo

! F' diagonalization now
! xnano will contain (X^-1)*C
! FFR: NO IT DOES NOT; it contains C....
         do i=1,M
            RMM(M15+i-1)=0.D0
            RMM(M13+i-1)=0.D0
         enddo

         do i=1,MM
            rmm5(i)=RMM(M5+i-1)
         enddo
         rmm15=0
         xnano=0


!        ESSL OPTION
#        ifdef  essl
         call DSPEV(1,RMM(M5),RMM(M13),X(1,M+1),M,M,RMM(M15),M2)
#        endif

!        LAPACK OPTION
#        ifdef pack
         call dspev('V','L',M,RMM5,RMM(M13),Xnano,M,RMM15,info)
#        endif

         do i =1,M
         do j=1,M
            X(i,M+j)=xnano(i,j)
         enddo
         enddo

! Recover C from (X^-1)*C
         do i=1,MM
!            RMM(M5+i-1)=rmm5(i)
         enddo

         do i=1,M
         do j=1,M
            X(i,M2+j)=0.D0
            do k=1,M
               X(i,M2+j)=X(i,M2+j)+X(i,k)*X(k,M+j)
            enddo
         enddo
         enddo

         call g2g_timer_stop('initial guess')


! Density Matrix
         kk=0
         do k=1,NCO
         do i=1,M
            kk=kk+1
            MO_coef_at(kk)=X(i,M2+k)
         enddo
         enddo

         kk=0
         do j=1,M
         do i=j,M
            kk=kk+1
            RMM(kk)=0.D0

            if (i.eq.j) then
!              one factor of 2 for alpha+beta
               ff=2.D0
            else
!              another factor of 2 for direct triangular sum (j>i)
!              w/real basis
               ff=4.D0
            endif
!
            do k=1,NCO
               RMM(kk)=RMM(kk)+ff*X(i,M2+k)*X(j,M2+k)
            enddo
         enddo
         enddo
!
         call g2g_timer_sum_stop('initial guess')
      endif
      deallocate(rmm5,rmm15)

!     End of Starting guess (No MO , AO known)-------------------------------
   END SUBROUTINE starting_guess_old



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   subroutine xPut( xMat, nCol, rMat )
   real*8 , intent(inout) :: xMat(:,:)
   integer, intent(in)    :: nCol
   real*8 , intent(in)    :: rMat(:,:)
   integer                :: ii, jj
 
   if ( size(xMat,1) /= size(rMat,1) ) then
      print*, "ERROR IN xPut - wrong number of rows"; stop
   end if

   if ( size(xMat,2) > (nCol + size(rMat,2)) ) then
      print*, "ERROR IN xPut - wrong number of cols"; stop
   end if
   
   do jj = nCol, size(rMat,2)
   do ii = 1, size(rMat,1)
      xMat(ii,jj) = rMat(ii,jj)
   end do
   end do

   end subroutine xPut

!------------------------------------------------------------------------------!
   subroutine xGet( xMat, nCol, rMat )
   real*8 , intent(in)    :: xMat(:,:)
   integer, intent(in)    :: nCol
   real*8 , intent(inout) :: rMat(:,:)
   integer                :: ii, jj

   if ( size(xMat,1) /= size(rMat,1) ) then
      print*, "ERROR IN xGet - wrong number of rows"; stop
   end if

   if ( size(xMat,2) < (nCol + size(rMat,2)) ) then
      print*, "ERROR IN xGet - wrong number of cols"; stop
   end if

   do jj = nCol, size(rMat,2)
   do ii = 1, size(rMat,1)
      rMat(ii,jj) = xMat(ii,jj)
   end do
   end do
   
   end subroutine xGet



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   SUBROUTINE COPY_VEC(VEC,DIM_VEC,POINTER_RMM)
!     subrutina temporal para empezar a romper RMM
!     copia el vector VEC a RMM posicion POINTER_RMM

      use garcha_mod, ONLY: RMM
      IMPLICIT NONE
      integer, intent(in) :: DIM_VEC,POINTER_RMM
      real*8, dimension(DIM_VEC), intent(in) :: VEC
      integer :: i

      do i=1, DIM_VEC
         RMM(POINTER_RMM+i-1)=VEC(i)
      end do

   END SUBROUTINE COPY_VEC



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
end module tmpaux_SCF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
