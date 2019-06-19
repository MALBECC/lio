!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% MAGNUS.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains the routine for an N-order magnus propagation.            !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module propagators
#include "complex_type.fh"
contains

subroutine predictor(F1a, F1b, FON, rho2, factorial, Xmat, Xtrans, timestep, &
                     time, M_in, MTB, dim3)
   ! This routine receives: F1a, F1b, rho2
   ! And gives: F5 = F(t+(deltat/2))
   use garcha_mod   , only: NBCH, rhoalpha, rhobeta, OPEN, r, d, natom,      &
                            ntatom, Iz, MEMO, Fmat_vec, Fmat_vec2, Ginv_vec, &
                            Hmat_vec, Gmat_vec, Pmat_vec
   use field_subs   , only: field_calc
   use mathsubs     , only: basechange
   use faint_cpu    , only: int3lu
   use fockbias_subs, only: fockbias_apply
   use basis_data   , only: M, Md
   use typedef_cumat, only: cumat_r, cumat_x

   implicit none
   ! MTB is only used in DFTB, it equals 0 otherwise.
   integer         , intent(in)    :: M_in, dim3, MTB
   TDCOMPLEX       , intent(in)    :: rho2(M_in,M_in,dim3)
   double precision, intent(in)    :: timestep, time, factorial(NBCH)
   type(cumat_r)   , intent(in)    :: Xmat
   type(cumat_x)   , intent(in)    :: Xtrans
   double precision, intent(inout) :: F1a(M_in,M_in,dim3), F1b(M_in,M_in,dim3),&
                                      FON(M_in,M_in,dim3)
   TDCOMPLEX, allocatable :: rho4(:,:,:), rho2t(:,:,:)
   integer :: i,j,k,kk, M2
   double precision :: E2, tdstep1, Ex, E1
   double precision, allocatable :: F3(:,:,:), FBA(:,:,:)

   allocate(rho4(M_in,M_in,dim3), rho2t(M_in,M_in,dim3), F3(M_in,M_in,dim3), &
            FBA(M_in,M_in,dim3))
   M2 = 2 * M

   ! Initializations and defaults
   ! tdstep of the predictor is 0.5 * tdstep_magnus
   ! F1a and F1b are used to extrapolate F3, then F3 is used to propagate rho.
   ! Afterwards, rho4 is copied into Pmat_vec in order to calculate F5.
   tdstep1 = timestep * 0.50D0
   F3      = 1.75D0 * F1b - 0.75D0 * F1a
   rho2t   = rho2

   call magnus(F3(:,:,1), rho2(:,:,1), rho4(:,:,1), M_in, NBCH, tdstep1, &
               factorial)
   if (OPEN) then
      call magnus(F3(:,:,2), rho2(:,:,2), rho4(:,:,2), M_in, NBCH, tdstep1, &
                  factorial)
      rho2t = rho4
      call Xtrans%change_base(rho2t(:,:,1), 'dir')
      call Xtrans%change_base(rho2t(:,:,2), 'dir')

      call sprepack_ctr('L', M, rhoalpha, rho2t(MTB+1:MTB+M,MTB+1:MTB+M,1))
      call sprepack_ctr('L', M, rhobeta , rho2t(MTB+1:MTB+M,MTB+1:MTB+M,2))
      Pmat_vec = rhoalpha + rhobeta
   else
      rho2t = rho4
      call Xtrans%change_base(rho2t(:,:,1), 'dir')
      call sprepack_ctr('L', M, Pmat_vec, rho2t(MTB+1:MTB+M,MTB+1:MTB+M,1))
   endif

   call int3lu(E2, Pmat_vec, Fmat_vec2, Fmat_vec, Gmat_vec, Ginv_vec, &
               Hmat_vec, open, MEMO)
   call g2g_solve_groups(0, Ex, 0)
   call field_calc(E1, time, Pmat_vec, Fmat_vec2, Fmat_vec, r, d, &
                   Iz, natom, ntatom, open)

   ! This is performed to recover TB terms from FON. If not, TB terms
   ! in FON become zero.
   FBA(:,:,1) = FON(:,:,1)
   call spunpack('L', M, Fmat_vec, FBA(MTB+1:MTB+M,MTB+1:MTB+M,1))
   call fockbias_apply(time, FBA(MTB+1:MTB+M,MTB+1:MTB+M,1))
   FON(:,:,1) = FBA(:,:,1)
   call Xmat%change_base(FON(:,:,1), 'dir')

   if (OPEN) then
      FBA(:,:,2) = FON(:,:,2)
      call spunpack('L', M, Fmat_vec2, FBA(MTB+1:MTB+M,MTB+1:MTB+M,2))
      call fockbias_apply(time, FBA(MTB+1:MTB+M,MTB+1:MTB+M,2))
      FON(:,:,2) = FBA(:,:,2)
      call Xmat%change_base(FON(:,:,2), 'dir')
   end if

   deallocate(rho4, rho2t, F3, FBA)
end subroutine predictor

subroutine magnus(Fock, RhoOld, RhoNew, M, N, dt, factorial)
   ! Propagator based in Baker-Campbell-Hausdorff Nth-order formula.
   ! Everything is calculated in the ortonormal basis.
   ! Input : Fock(t+(deltat/2)), rho(t)
   ! Output: rho6 = rho(t+deltat)
   use typedef_cumat, only: cumat_x

   implicit none
   integer     , intent(in)  :: M, N
   real(kind=8), intent(in)  :: Fock(M,M), factorial(N), dt
   TDCOMPLEX   , intent(in)  :: RhoOld(M,M)
   TDCOMPLEX   , intent(out) :: RhoNew(M,M)

   integer                :: icount, jcount
   type(cumat_x)          :: omega, comm_prev, comm_next, rho_m
   TDCOMPLEX              :: alpha, beta
   TDCOMPLEX, allocatable :: Omega1(:,:)
   TDCOMPLEX, parameter   :: ICMPLX = CMPLX(0.0D0, 1.0D0)


   allocate(Omega1(M,M))
   do icount = 1, M
   do jcount = 1, M
      Omega1(icount,jcount) = - ICMPLX * Fock(icount,jcount) * dt
   enddo
   enddo

   call omega%init(M, Omega1, .true.)
   call rho_m%init(M, rhoOld, .true.)
   call comm_prev%init(M, rhoOld, .true.)
   call comm_next%init(M, rhoOld, .true.)

   ! Density matrix propagation
   alpha  = (1.0D0, 0.0D0)
   do icount = 1, N
      alpha = CMPLX(factorial(icount), 0.0D0)
      beta  = CMPLX(0.0D0, 0.0D0)
      call comm_next%mat_mul(omega    , comm_prev, alpha, beta)

      beta  = CMPLX(-1.0D0, 0.0D0)
      call comm_next%mat_mul(comm_prev, omega    , alpha, beta)
      
      beta  = CMPLX(1.0D0, 0.0D0)
      call rho_m%add_mat(comm_next, beta)

      call comm_next%exchange(comm_prev)
   enddo

   ! Stores the new density.
   call rho_m%get(rhoNew)

   call omega%destroy()
   call rho_m%destroy()
   call comm_next%destroy()
   call comm_prev%destroy()
   deallocate(Omega1)
end subroutine magnus

end module propagators
