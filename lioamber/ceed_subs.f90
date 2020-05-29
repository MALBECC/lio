#include "datatypes/datatypes.fh"
module ceed_subs
   implicit none
contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ceed_init(M, open_shell, r, d, natom, ntatom, propagator)
! This subroutine initialize the variables for CEED calculations
   use ceed_data, only: ceed_calc, dip_ceed_op, d2ip_ceed_op,            &
                         fock_ceed_op, d2dip_ceed, Xmat_ceed, k_ceed
   use faint_cpu , only: intfld

   implicit none
   logical, intent(in)  :: open_shell
   integer, intent(in)  :: M
   integer, intent(in)  :: propagator
   integer, intent(in)  :: natom
   integer, intent(in)  :: ntatom
   LIODBLE, intent(in)  :: r(ntatom,3)
   LIODBLE, intent(in)  :: d(natom,natom)
   integer              :: ii
   integer              :: MM
   LIODBLE              :: vec_aux(3)
   LIODBLE, allocatable :: dip_array(:)
   LIODBLE, allocatable :: dip_mat_aux(:,:)

   if (.not. ceed_calc) return

   if (propagator /= 1) then
      write(*,*) "CEED calculation can only be performed with Verlet propagator"
      write(*,*) "Stopping LIO"
      stop
   end if

   MM = M*(M+1)/2

   allocate(dip_ceed_op(3),dip_array(MM), dip_mat_aux(M, M))

   if (.not.open_shell) then
      allocate(fock_ceed_op(1), d2ip_ceed_op(3,1), d2dip_ceed(3,1))
   else
      allocate(fock_ceed_op(2), d2ip_ceed_op(3,2), d2dip_ceed(3,2))
   end if

   do ii=1, 3
      dip_array    = 0.0d0
      dip_mat_aux  = 0.0d0
      vec_aux      = 0.0d0
      vec_aux(ii)  = 1.0d0
      call intfld(dip_array, dip_array, r, d, natom, ntatom, .false., 1.0d0,   &
                  vec_aux(1), vec_aux(2), vec_aux(3))
      call spunpack('L', M, dip_array, dip_mat_aux)
      call dip_ceed_op(ii)%Sets_data_AO(dip_mat_aux)
      call dip_ceed_op(ii)%BChange_AOtoON(Xmat_ceed, M)
   end do

end subroutine ceed_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ceed_fock_calculation(fock, rho_aux, M, t_step, dim3, open_shell)
!This subroutine includes the CEED term in the evolution of the density matrix.
!At the end of the subroutine rho_aux store the commutator -i[H_CEED,\rho],
!were:
!              H_CEED =  -i A_ceed d^2<\mu>/dt^2 . [\mu , H]
!d^2<\mu>/dt^2 is approximated as:
!               d^2<\mu>/dt^2 = - Tr(\rho . [[\mu,H],H])
!and:
!             A_ceed = 2/3 \alpha/c^2 = 2.592668788247536d-7

   use ceed_data, only: ceed_calc, dip_ceed_op, d2ip_ceed_op, d2dip_ceed,     &
                         A_ceed, fock_ceed_op, k_ceed,ceed_td_step
   implicit none

   logical   , intent(in)    :: open_shell
   integer   , intent(in)    :: M
   integer   , intent(in)    :: dim3
   integer   , intent(in)    :: t_step
   LIODBLE   , intent(in)    :: fock(M,M,dim3)
   TDCOMPLEX , intent(inout) :: rho_aux(M,M,dim3)
   integer                :: ii, jj, kk, ss
   LIODBLE                :: aux1
   LIODBLE                :: aux_mat1(M,M,dim3)
   LIODBLE                :: aux_mat2(M,M,dim3)
   LIODBLE                :: aux_mat3(M,M,dim3)
   TDCOMPLEX              :: aux_mat4(M,M,dim3)

   if (.not.ceed_calc) return

   aux_mat3 = 0.0d0
   aux_mat4 = 0.0d0

   if (t_step > ceed_td_step) then
      do ii = 1,3
      do ss = 1,dim3
         call dip_ceed_op(ii)%Commut_data_r(fock(:,:,ss), aux_mat1(:,:,ss), M)
         call d2ip_ceed_op(ii,ss)%Sets_data_ON(aux_mat1(:,:,ss))
         call d2ip_ceed_op(ii,ss)%Commut_data_r(fock(:,:,ss),aux_mat2(:,:,ss),M)
         d2dip_ceed(ii,ss) = 0.0d0
         do jj=1, M
         do kk=1, M
            aux1 = dble(rho_aux(jj,kk,ss)*aux_mat2(kk,jj,ss))
            d2dip_ceed(ii,ss) = d2dip_ceed(ii,ss) - aux1
         end do
         end do

      end do
         aux1     = d2dip_ceed(ii,1)
         if (open_shell) aux1 = aux1 + d2dip_ceed(ii,2)
         aux_mat3 = aux_mat3 + k_ceed * A_ceed * aux1 * aux_mat1
      end do

      do ss = 1,dim3
         call fock_ceed_op(ss)%Sets_data_ON(aux_mat3(:,:,ss))
         call fock_ceed_op(ss)%Commut_data_c(rho_aux(:,:,ss),                  &
                                             aux_mat4(:,:,ss), M)
      end do

   end if

   rho_aux = aux_mat4


end subroutine ceed_fock_calculation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

end module ceed_subs
