!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine simple_guess( atom_pos, dens_mat )
!
! Is a subroutine like this really necessary? Couldn't it be replaced by
! enabling a 0 pass on the SCF cycle where only int1 is computed?
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   
   use garcha_mod, only: M, RMM, VCINP, primera, X, Md, NCO
   implicit none
   real*8, intent(out) :: mocoef_at(M,M)

   integer :: info
   real*8, intent(out) ::xnano
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

      write(666,*) "RMM en M11"
      do ii=1,M
         kk = (M2-ii)*(ii-1)/2
         do jj=ii+1,M
            write(666,*) RMM(M11+jj+kk-1)
         enddo
      enddo
      write(666,*) 

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

         write(666,*) " INFO = ", info

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


end subroutine simple_guess
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
