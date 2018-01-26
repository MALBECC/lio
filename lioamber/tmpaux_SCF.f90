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
!  obtain_new_P: improves fock guess with the convergence accelerator and
!                diagonalizes it to obtains the coeficientes. This two tasks
!                should be separated.
!
!  density: performs the construction of density matrix from the orbital 
!           coefficients.
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
   SUBROUTINE starting_guess(xnano)
      use garcha_mod, ONLY: RMM, VCINP, primera, M, X, Md, NCO, MO_coef_at
      IMPLICIT NONE
      integer :: info
      real*8, dimension (M,M), intent(inout)::xnano
      real*8, dimension (:), ALLOCATABLE :: rmm5,rmm15
      integer :: M1,M2,M3, M5, M7, M9, M11, M13, M15, M17, M18, MM, MMd 
               ! temporales hasta q rompamos RMM
      integer :: i,j,k,kk !auxiliares
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
            RMM(M5+i-1)=rmm5(i)
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
   END SUBROUTINE starting_guess



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!carlos: sacando fock de la rutina
#ifdef CUBLAS
   subroutine obtain_new_P( M_in, NCO_in, niter, damp, good, fock, rho, &
                          & morb_energy, morb_coefat, devPtrX, devPtrY)
      use cublasmath, only : cumxp_r, cumfx, cumxtf, cu_calc_fock_commuts
#else
   subroutine obtain_new_P( M_in, NCO_in, niter, damp, good, fock, rho, &
                          & morb_energy, morb_coefat, Xmat, Ymat )
      use mathsubs  , only : basechange_gemm
#endif

      use garcha_mod    , only: DIIS, NCO, ndiis, hybrid_converg, good_cut
      use linear_algebra, only: matrix_diagon
      use converger_subs, only: converger_init, conver

      implicit none
      integer  , intent(in)    :: M_in
      integer  , intent(in)    :: NCO_in
      integer  , intent(in)    :: niter
      real*8   , intent(in)    :: damp
      real*8   , intent(in)    :: good
      real*8   , intent(inout)    :: fock(M_in, M_in)
      real*8   , intent(inout) :: rho(M_in,M_in)
      real*8   , intent(inout) :: morb_energy(M_in)
      real*8   , intent(inout) :: morb_coefat(M_in, M_in)
#     ifdef  CUBLAS
      integer*8, intent(in)    :: devPtrX
      integer*8, intent(in)    :: devPtrY
#     else
      real*8   , intent(in)    :: Xmat(M_in,M_in)
      real*8   , intent(in)    :: Ymat(M_in,M_in)
#     endif
      real*8, allocatable :: eigen_vecs(:,:)
      integer :: ii, jj


!  Convergence acceleration 
!------------------------------------------------------------------------------!
      call g2g_timer_sum_start('SCF acceleration')
      if (niter==1) call converger_init( M_in, ndiis, damp, DIIS, hybrid_converg )

#     ifdef CUBLAS
         call conver(niter, good, good_cut, M_in, rho, fock, devPtrX, devPtrY )
#     else
         call conver(niter, good, good_cut, M_in, rho, fock, Xmat, Ymat)
#     endif
      call g2g_timer_sum_pause('SCF acceleration')


!  Fock(ON) diagonalization, base change of coeficients ( (X^-1)*C ) and
!  construction of new density matrix
!------------------------------------------------------------------------------!
      if ( allocated(eigen_vecs) ) deallocate(eigen_vecs)
      allocate( eigen_vecs(M_in,M_in) )


      call g2g_timer_start('SCF - Fock Diagonalization')
      call g2g_timer_sum_start('SCF - Fock Diagonalization (sum)')
      call matrix_diagon( fock, eigen_vecs, morb_energy )
      call g2g_timer_sum_pause('SCF - Fock Diagonalization (sum)')
      call g2g_timer_stop('SCF - Fock Diagonalization')

      call g2g_timer_start('SCF - MOC base change')
      call g2g_timer_sum_start('SCF - MOC base change (sum)')
#     ifdef CUBLAS
         call cumxp_r( eigen_vecs, devPtrX, morb_coefat, M_in)
#     else
         morb_coefat = matmul( Xmat, eigen_vecs )
#     endif
      call g2g_timer_sum_pause('SCF - MOC base change (sum)')
      call g2g_timer_stop('SCF - MOC base change')

      if ( allocated(eigen_vecs) ) deallocate(eigen_vecs)

      call density( M_in, NCO_in, morb_coefat, rho )

end subroutine obtain_new_P



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine density(M, NCO, X, C)
    use garcha_mod, only : OPEN   
    ! M:    number of atomic basis functions.
    ! NCO:  number of occupied orbitals.
    ! X:    matrix containing Fock coefficients.
    ! C:    matrix containing output density matrix.
    implicit none
    integer,intent(in)  :: M, NCO
    real*8, intent(in)  :: X(M, M)
    real*8, intent(out) :: C(M, M)

    real*8, allocatable :: A(:, :)
    real*8              :: factor
    integer             :: DOSM, i, j
#ifdef CUBLAS
    integer*8 devPtrA, devPtrC
    integer   sizeof_real
    parameter(sizeof_real=8)
    integer   stat
    external  CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
    external  CUBLAS_SHUTDOWN, CUBLAS_ALLOC,CUBLAS_FREE
    integer   CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
    integer   CUBLAS_INIT
#endif

    ! This factor take into account occupation for closed-shell
    ! systems.
    factor = 2.0D0 
    if (OPEN) factor = 1.0D0

    allocate(A(M,NCO))
    do i=1, M
        do j=1, NCO
             A(i, j) = X(i, j)
        enddo
    enddo
    C = 0.0D0

#ifdef CUBLAS
    ! Allocates matrix A on the device.
    stat = CUBLAS_ALLOC(M*NCO, sizeof_real, devPtrA)
    if (stat.NE.0) then
        write(*,*) "density: A-Matrix allocation failed."
        call CUBLAS_SHUTDOWN
        stop
    endif

    ! Allocate matrix C on the device
    stat = CUBLAS_ALLOC(M*M, sizeof_real, devPtrC)
    if (stat.NE.0) then
        write(*,*) "density: C-Matrix allocation failed."
        call CUBLAS_SHUTDOWN
        stop
    endif
 
    ! Copy A in the device
    stat = CUBLAS_SET_MATRIX(M, NCO, sizeof_real, A, M, devPtrA, M)
    if (stat.NE.0) then
        write(*,*) "matrix copy to the device failed -density"
        call CUBLAS_SHUTDOWN
        stop 
    endif

    ! Peforms the matrix multiplication, obtaining C as A*A.
    call CUBLAS_DGEMM('N', 'T', M, M, NCO, factor, devPtrA, M, devPtrA, M, &
                      0.0D0, devPtrC, M)

    ! Copies C to host.
    stat = CUBLAS_GET_MATRIX(M, M, sizeof_real, devPtrC, M, C, M)

    ! Finalizes variables.
    call CUBLAS_FREE(devPtrA)
    call CUBLAS_FREE(devPtrC)

#else
    ! Obtains C as A*A.
    call DGEMM('N', 'T', M, M, NCO, factor, A, M, A, M, 0.0D0, C, M)
#endif
    
    ! Multiplies by 2 all non diagonal elements.
    do i=1,M
        do j=1,i-1
            C(i,j) = 2.0D0 * C(i,j)
        enddo
        do j=i+1,M
            C(i,j) = 2.0D0 * C(i,j)
        enddo
    enddo

    deallocate(A)
    return
end subroutine density


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
