!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!     CASE OF NO STARTING GUESS PROVIDED, 1 E FOCK MATRIX USED
!     FCe = SCe; (X^T)SX = 1
!     F' = (X^T)FX
!     => (X^-1*C)^-1 * F' * (X^-1*C) = e
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine starting_guess( Nmat, Nvec, hmat_vec, Xmat, densat_vec )
   use garcha_mod, ONLY:  M, X, Md, NCO, MO_coef_at
   use liosubs_math, only: transform
   implicit none
   integer, intent(in)    :: Nmat
   integer, intent(in)    :: Nvec
   real*8 , intent(in)    :: Xmat(Nmat,Nmat)
   real*8 , intent(inout) :: hmat_vec(Nvec)
   real*8 , intent(inout) :: densat_vec(Nvec)
   
   real*8, allocatable   :: morb_coefat(:,:)
   real*8, allocatable   :: morb_coefon(:,:)
   real*8, allocatable   :: rmm5(:), rmm15(:)
   real*8, allocatable   :: hmat(:,:), fockon(:,:), fockon_vec(:)
   real*8, allocatable   :: eigval_vec(:)
   integer :: ii,jj,i,j,k,kk
   integer :: M1, M2, info
   real*8  :: ff

   call g2g_timer_start('initial guess')
   call g2g_timer_sum_start('initial guess')


   allocate( rmm5(Nvec), rmm15(Nvec), fockon_vec(Nvec) )
   allocate( eigval_vec(Nvec) )
   allocate( morb_coefat(Nmat, Nmat), morb_coefon(Nmat, Nmat) )

!  Calculate F' in RMM(M5)
   allocate( hmat(Nmat,Nmat) )
   allocate( fockon(Nmat, Nmat) )
   call spunpack('L', Nmat, hmat_vec, hmat )
   fockon = transform( hmat, Xmat )
   call sprepack('L', Nmat, fockon_vec, fockon )
 

! F' diagonalization now
! xnano will contain (X^-1)*C
! FFR: NO IT DOES NOT; it contains C...
   do i=1,M
      rmm15(i)=0.D0
      eigval_vec(i) = 0.0d0
   enddo

   do i=1,Nvec
      rmm5(i)=fockon_vec(i)
   enddo

   morb_coefon(:,:) = 0.0d0

!  Inputs and outputs are the same
   call dspev( 'V', 'L', M, RMM5, eigval_vec, morb_coefon, M, RMM15, info )


! Recover C from (X^-1)*C
   do i=1,M
   do j=1,M
      morb_coefat(i,j) = 0.0d0
      do k=1,M
         morb_coefat(i,j) = morb_coefat(i,j) + X(i,k) * morb_coefon(k,j)
      enddo
   enddo
   enddo

! Density Matrix
   kk=0
   do k=1,NCO
   do i=1,M
      kk=kk+1
      MO_coef_at(kk)=morb_coefat(i,k)
   enddo
   enddo

   
   kk=0
   do j=1,M
   do i=j,M
   
      kk=kk+1
      densat_vec(kk)=0.D0
      
      if (i.eq.j) then
!        one factor of 2 for alpha+beta
         ff=2.D0
      else
!        another factor of 2 for direct triangular sum (j>i) w/real basis
         ff=4.D0
      endif
!
      do k=1,NCO
         densat_vec(kk)=densat_vec(kk)+ff*morb_coefat(i,k)*morb_coefat(j,k)
      enddo
      
   enddo
   enddo
 
   M1 = Nmat
   M2 = 2*Nmat
   do jj = 1, M
   do ii = 1, M
      X(ii,M1+jj)       = morb_coefon(ii,jj)
      X(ii,M2+jj)       = morb_coefat(ii,jj)
   enddo
   enddo


   deallocate(rmm5,rmm15)

   call g2g_timer_stop('initial guess')
   call g2g_timer_sum_stop('initial guess')
end subroutine starting_guess
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
