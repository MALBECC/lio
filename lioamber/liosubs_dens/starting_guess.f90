!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!     CASE OF NO STARTING GUESS PROVIDED, 1 E FOCK MATRIX USED
!     FCe = SCe; (X^T)SX = 1
!     F' = (X^T)FX
!     => (X^-1*C)^-1 * F' * (X^-1*C) = e
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine starting_guess( Nmat, Nvec, hmat_vec, fockat_vec, eigval_vec,       &
                         & eigvec_mat, densat_vec )
   use garcha_mod, ONLY: RMM, M, X, Md, NCO, MO_coef_at
   implicit none
   integer, intent(in)    :: Nmat
   integer, intent(in)    :: Nvec
   real*8 , intent(inout) :: hmat_vec(Nvec)
   real*8 , intent(inout) :: fockat_vec(Nvec)
   real*8 , intent(inout) :: eigval_vec(Nvec)
   real*8 , intent(inout) :: eigvec_mat(Nmat,Nmat)
   real*8 , intent(inout) :: densat_vec(Nvec)


   real*8, allocatable   :: rmm5(:), rmm15(:)
   integer :: MM, MMd, M1, M2, M3, M5, M7, M9, M11, M13, M15
   integer :: ii,jj,i,j,k,kk
   integer :: info
   real*8  :: ff

   call g2g_timer_start('initial guess')
   call g2g_timer_sum_start('initial guess')


   MM  = M*(M+1)/2
   MMd = Md*(Md+1)/2

   M1  = 1         ! first P
   M2  = 2*M
   M3  = M1  + MM  ! now Pnew
   M5  = M3  + MM  ! now S, F also uses the same position after S was used
   M7  = M5  + MM  ! now G
   M9  = M7  + MMd ! now Gm
   M11 = M9  + MMd ! now H
   M13 = M11 + MM  ! W ( eigenvalues ), also this space is used in least squares
   M15 = M13 + M   ! aux ( vector for ESSl)

   allocate(rmm5(MM),rmm15(mm))


!  Calculate F' in RMM(M5)

   do i=1,M
   do j=1,M
      X(i,M+j)=0.D0
      do k=1,j
         X(i,M+j)=X(i,M+j)+X(k,i)*hmat_vec(j+(M2-k)*(k-1)/2)
      enddo
      do k=j+1,M
         X(i,M+j)=X(i,M+j)+X(k,i)*hmat_vec(k+(M2-j)*(j-1)/2)
      enddo
   enddo
   enddo

   kk=0
   do j=1,M
   do i=j,M
      kk=kk+1
      fockat_vec(kk)=0.D0
      do k=1,j
         fockat_vec(kk)=fockat_vec(kk)+X(i,M+k)*X(k,j)
      enddo
   enddo
   enddo

! F' diagonalization now
! xnano will contain (X^-1)*C
   do i=1,M
      RMM(M15+i-1)=0.D0
      eigval_vec(i) = 0.0d0
   enddo

   do i=1,MM
      rmm5(i)=fockat_vec(i)
   enddo
   rmm15=0
   eigvec_mat=0
   
   call dspev('V','L',M,RMM5,eigval_vec,eigvec_mat,M,RMM15,info)
   write(666,*) " INFO = ", info

   do i =1,M
   do j=1,M
      X(i,M+j)=eigvec_mat(i,j)
   enddo
   enddo

! Recover C from (X^-1)*C
   do i=1,M
   do j=1,M
      X(i,M2+j)=0.D0
      do k=1,M
         X(i,M2+j)=X(i,M2+j)+X(i,k)*X(k,M+j)
      enddo
   enddo
   enddo

   
   write(666,*) "Coeficientes en base atomica"
   do ii=1,M
   do jj=1,M
      write(666,*) ii, jj, X(ii,M2+jj)
   enddo
   enddo
   write(666,*)


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
         densat_vec(kk)=densat_vec(kk)+ff*X(i,M2+k)*X(j,M2+k)
      enddo
      
   enddo
   enddo
   
   deallocate(rmm5,rmm15)

   call g2g_timer_stop('initial guess')
   call g2g_timer_sum_stop('initial guess')
end subroutine starting_guess
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
