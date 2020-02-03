subroutine PCG_solve(bvec,Coef,E,X,M,NCO,Nvirt,Ndim)
use excited_data, only: fittExcited
use garcha_mod,   only: PBE0
   implicit none

   integer, intent(in) :: M, NCO, Nvirt, Ndim
   double precision, intent(in) :: bvec(Ndim), E(M), Coef(M,M)
   double precision, intent(inout) :: X(Ndim)

   integer :: iter, maxIter
   logical :: conv
   double precision :: beta, alpha
   double precision, dimension(:), allocatable :: R, Pk, Mprec, ApIA
   double precision, dimension(:,:), allocatable :: Pmat, F2e, Fxc
   double precision, dimension(:,:), allocatable :: CopyP

!  START PRECONDITINED CONJUGATE GRADIENT
   maxIter = 50; conv = .false.

!  INITIAL GUESS: Xo = 0
   allocate(R(Ndim))
   R = bvec

!  CALCULATE PRECONDITIONED M^(-1)
   allocate(Mprec(Ndim))
   call Prec_calculate(E,Mprec,M,NCO,Ndim)

   allocate(Pk(Ndim)); Pk = 0.0D0
   call Pbeta_calc(R,Mprec,beta,Pk,Ndim)

   allocate(Pmat(M,M),F2e(M,M),Fxc(M,M))
   allocate(ApIA(Ndim),CopyP(M,M))
   X = 0.0D0

   write(*,"(1X,A)") "Start PCG loop"
   do iter=1,maxIter

!     CONVERT TRIAL VECTORS TO AO BASIS
      call VecToMat(Pk,Pmat,Coef,Ndim,NCO,M)

!     CALCULATE TWO ELECTRON PART
      if ( .not. fittExcited ) then
         call g2g_timer_start("Fock 2e LR")
         call g2g_calculate2e(Pmat,F2e,1)
         F2e = (F2e+transpose(F2e))
         call g2g_timer_stop("Fock 2e LR")
      elseif ( fittExcited .and. (.not. PBE0) ) then
         call g2g_timer_start("Fock 2e LR")
         call calc2eFITT(Pmat,F2e,M)
         call g2g_timer_stop("Fock 2e LR")
      else
         print*, "Error in 2 Electron Repulsion Integrals"
         print*, "Check PBE0 and fittExcited"
         stop
      endif

!     CALCULATE XC PART
      Fxc = 0.0d0
      call g2g_calculateg(Pmat,Fxc,2)

!     OBTAIN FOCK TOTAL AND ADD TERM (Ea-Ei)Pk AND
!     CHANGE BASIS AO -> MO
      call Ap_calculate(F2e,Fxc,Pk,E,ApIA,M,NCO,Nvirt,Ndim)
     
!     CALCULATE ALPHA
      call Alpha_calc(Pk,ApIA,alpha,Ndim)
     
!     NEW X VECTOR
      X = X + alpha * Pk
     
!     NEW R VECTOR
      R = R - alpha * APIA
     
!     CHECK CONVERGENCE
      call error(R,conv,Ndim,iter)
      if (conv .eqv. .true.) then
         write(*,"(1X,A,1X,I2,1X,A)") "Convergence achieved in",iter,"itertaions"
         exit
      endif

!     GET NEW BETA AND Pk
      call Pbeta_calc(R,Mprec,beta,Pk,Ndim)

   enddo ! ENDDO LOOP PCG

   deallocate(R,Mprec,Pk,Pmat,F2e,Fxc,ApIA,CopyP)
end subroutine PCG_solve

subroutine Prec_calculate(Ener,Mprec,M,NCO,Ndim)
   implicit none

   integer, intent(in) :: M, NCO, Ndim
   double precision, intent(in) :: Ener(M)
   double precision, intent(out) :: Mprec(Ndim)

   integer :: ii, jj, Nvirt, NCOc, pos
   double precision :: temp

   Nvirt = M - NCO
   NCOc = NCO + 1

   do ii=1,NCO
   do jj=1,Nvirt
      pos = (ii-1) * Nvirt + jj
      temp = Ener(jj+NCO) - Ener(NCOc-ii)
      temp = 1.0D0 / temp
      Mprec(pos) = temp
   enddo
   enddo
end subroutine Prec_calculate

subroutine Pbeta_calc(R,M,beta,P,N)
   implicit none

   integer, intent(in) :: N
   double precision, intent(in) :: R(N), M(N)
   double precision, intent(inout) :: P(N)
   double precision, intent(out) :: beta

   integer :: ii
   double precision :: temp

   temp = 0.0D0
   do ii=1,N
      temp = temp + R(ii) * R(ii) * M(ii)
   enddo
   beta = 1.0D0 / temp

   do ii=1,N
     P(ii) = P(ii) + beta * M(ii) * R(ii)
   enddo
end subroutine Pbeta_calc

subroutine VecToMat(Vec,Mat,Coef,Ndim,NCO,M)
use excited_data, only: Coef_trans

   implicit none

   integer, intent(in) :: Ndim, NCO, M
   double precision, intent(in) :: Vec(Ndim), Coef(M,M)
   double precision, intent(out) :: Mat(M,M)

   integer :: row, col, NCOc, Nvirt, pos
   double precision, dimension(:,:), allocatable :: scr

   Nvirt = M - NCO
   NCOc = NCO + 1

   Mat = 0.0D0
   do row=1,NCO
   do col=1,Nvirt
      pos = (row-1) * Nvirt + col
      Mat(NCOc-row,NCO+col) = Vec(pos)
   enddo
   enddo

   allocate(scr(M,M))
   call dgemm('N','N',M,M,M,1.0d0,Coef,M,Mat,M,0.0d0,scr,M)
   call dgemm('N','N',M,M,M,1.0d0,scr,M,Coef_trans,M,0.0d0,Mat,M)
   deallocate(scr)
end subroutine VecToMat

subroutine Ap_calculate(Fe,Fx,P,E,Ap,M,NCO,Nvirt,Ndim)
use excited_data, only: Cocc_trans, Cvir
   implicit none

   integer, intent(in) :: M, NCO, Nvirt, Ndim
   double precision, intent(in) :: Fe(M,M), E(M), P(Ndim)
   double precision, intent(inout) :: Fx(M,M)
   double precision, intent(out) :: Ap(Ndim)

   integer :: ii, jj, NCOc, pos
   double precision, allocatable :: Fp(:,:), scrA(:,:), scrB(:,:)

   ! Obtain total fock
   allocate(Fp(M,M))
   do ii=1,M
   do jj=ii,M
      Fx(jj,ii) = Fx(ii,jj)
   enddo
   enddo
   Fp = Fe + 2.0d0 * Fx

   ! Basis Change
   allocate(scrA(M,Nvirt),scrB(NCO,Nvirt))
   call dgemm('N','N',M,Nvirt,M,1.0d0,Fp,M,Cvir,M,0.0d0,scrA,M)
   call dgemm('N','N',NCO,Nvirt,M,1.0d0,Cocc_trans,NCO,scrA,&
              M,0.0d0,scrB,NCO)

!  Obtain A*p in MO Basis
   NCOc = NCO + 1
   do ii=1,NCO
   do jj=1,Nvirt
      pos = (ii-1) * Nvirt + jj
      Ap(pos) = scrB(NCOc-ii,jj) + ( E(NCO+jj) - E(NCOc-ii) )*P(pos)
   enddo
   enddo

   deallocate(scrA,scrB,Fp)
end subroutine Ap_calculate

subroutine Alpha_calc(P,A,alpha,N)
   implicit none

   integer, intent(in) :: N
   double precision, intent(in) :: P(N), A(N)
   double precision, intent(out) :: alpha

   integer :: ii
   double precision :: temp

   temp = 0.0D0
   do ii=1,N
     temp = temp + P(ii) * A(ii)
   enddo
   alpha = 1.0D0 / temp
end subroutine Alpha_calc

subroutine error(V,convergence,N,iter)
  implicit none

  integer, intent(in) :: N, iter
  double precision, intent(in) :: V(N)
  logical, intent(inout) :: convergence

  integer :: ii
  real*8 :: temp, tol

  tol = 1.0D-12
  temp = 0.0D0
  do ii=1,N
    temp = temp + V(ii)*V(ii)
  enddo

  write(*,8070) iter, temp, tol
  if( temp < tol ) convergence = .true.

8070  FORMAT(1X,"Iteration = ", I2,1X,&
      "Error (crit) =",ES9.2," (",ES9.2,")")
end subroutine error
