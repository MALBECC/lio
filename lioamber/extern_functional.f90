#include "datatypes/datatypes.fh"
! This module enables the usage of functionals from libxc

! VARIABLES:
! extern_functional = bool
! functional_id = identifier of functional in libxc ( see main page of libxc )
! HF = int array: 0 = false , 1 = true 
!      0 = Exact Hartree-Fock term ( FULL RANGE )
!      1 = Short-Range Hartree-Fock term ( SHORT RANGE )
!      2 = Long-Range Hartree-Fock term ( LONG RANGE )
! HF_fac = double arrays 
!      0 = factor of exact exchange
!      1 = factor of short range exchange
!      2 = factor of long range exchange
! screen = facotr short or long range, valid when HF[1] or HF[2] is true
!          if both cases are true, they use the same screen factor

module extern_functional_data
   implicit none

   ! Input Variables
   logical :: extern_functional = .false.
   integer :: functional_id     = -1

   ! Internal Variables
   integer :: HF(3) = (/0,0,0/)
   LIODBLE :: HF_fac(3) = (/0.0d0,0.0d0,0.0d0/)
   LIODBLE :: screen = 0.0d0
   logical :: libint_inited = .false.

   ! All HF energies components
   LIODBLE, dimension(:,:,:), allocatable :: FockHF_a0
   LIODBLE, dimension(:,:,:), allocatable :: FockHF_b0

end module extern_functional_data

module extern_functional_subs
   implicit none
contains
subroutine libint_init(c_raw,libint_recalc)
use fstsh_data, only: FSTSH
use excited_data, only: lresp
use extern_functional_data, only: HF, libint_inited
   implicit none

   LIODBLE, intent(in) :: c_raw(:,:) 
   integer, intent(in) :: libint_recalc

   integer :: ii
   logical :: need_libint, final_decision

   need_libint = .false.
   do ii=1,3
      if ( HF(ii)==1 ) need_libint = .true.
   enddo

   final_decision = need_libint .or. libint_inited
   final_decision = final_decision .or. lresp
   final_decision = final_decision .or. FSTSH

   if ( final_decision ) then
      call g2g_timer_sum_start('Libint init')
      call g2g_libint_init(c_raw,libint_recalc)
      call g2g_timer_sum_stop('Libint init')
      libint_inited = .true.
   endif
end subroutine libint_init

subroutine exact_exchange(rho_a0,rho_b0,fock_a0,fock_b0,M)
use garcha_mod, only: OPEN
use extern_functional_data, only: HF, HF_fac, FockHF_a0, FockHF_b0
   implicit none

   integer, intent(in) :: M
   LIODBLE, intent(in) :: rho_a0(M,M), rho_b0(M,M)
   LIODBLE, intent(inout) :: fock_a0(M,M), fock_b0(M,M)

   if ( allocated(FockHF_a0) ) deallocate(FockHF_a0)
   allocate(FockHF_a0(3,M,M)); FockHF_a0 = 0.0d0

   call g2g_timer_sum_start('All Exact Exchange Fock')
   if ( OPEN ) then 
      if ( allocated(FockHF_b0) ) deallocate(FockHF_b0)
      allocate(FockHF_b0(3,M,M)); FockHF_b0 = 0.0d0

      ! Exact Exchange Hartree Fock
      if ( HF(1)==1 ) then
         call g2g_exact_exchange_open(rho_a0,rho_b0, &
                        FockHF_a0(1,:,:),FockHF_b0(1,:,:),1)
         fock_a0 = fock_a0 - HF_fac(1) * FockHF_a0(1,:,:)
         fock_b0 = fock_b0 - HF_fac(1) * FockHF_b0(1,:,:)
      endif

      ! Short Range Exchange Hartree Fock
      if ( HF(2)==1 ) then
         call g2g_exact_exchange_open(rho_a0,rho_b0,&
                        FockHF_a0(2,:,:),FockHF_b0(2,:,:),2)
         fock_a0 = fock_a0 - HF_fac(2) * FockHF_a0(2,:,:)
         fock_b0 = fock_b0 - HF_fac(2) * FockHF_b0(2,:,:)
      endif

      ! Long Range Exchange Hartree Fock
      if ( HF(3)==1 ) then
         call g2g_exact_exchange_open(rho_a0,rho_b0,&
                        FockHF_a0(3,:,:),FockHF_b0(3,:,:),3)
         fock_a0 = fock_a0 - HF_fac(3) * FockHF_a0(3,:,:)
         fock_b0 = fock_b0 - HF_fac(3) * FockHF_b0(3,:,:)
      endif

  else ! CLOSED SHELL
      
      ! Exact Exchange Hartree Fock
      if ( HF(1)==1 ) then
         call g2g_exact_exchange(rho_a0,FockHF_a0(1,:,:),1)
         fock_a0 = fock_a0 - HF_fac(1) * FockHF_a0(1,:,:)
      endif

      ! Short Range Exchange Hartree Fock
      if ( HF(2)==1 ) then
         call g2g_exact_exchange(rho_a0,FockHF_a0(2,:,:),2)
         fock_a0 = fock_a0 - HF_fac(2) * FockHF_a0(2,:,:)
      endif

      ! Long Range Exchange Hartree Fock
      if ( HF(3)==1 ) then
         call g2g_exact_exchange(rho_a0,FockHF_a0(3,:,:),3)
         fock_a0 = fock_a0 - HF_fac(3) * FockHF_a0(3,:,:)
      endif
   endif ! END CLOSED SHELL
   call g2g_timer_sum_pause('All Exact Exchange Fock')
end subroutine exact_exchange

subroutine exact_energies(rho_a0,rho_b0,E1,E2,E3,M)
use garcha_mod, only: OPEN
use extern_functional_data, only: FockHF_a0, FockHF_b0, HF_fac
   implicit none

   integer, intent(in) :: M
   LIODBLE, intent(out) :: E1, E2, E3
   LIODBLE, intent(in) :: rho_a0(M,M), rho_b0(M,M)

   integer :: ii, jj

   call g2g_timer_sum_start("All Exact Exchange Energy")
   E1 = 0.0d0; E2 = 0.0d0; E3 = 0.0d0
   do ii=1,M
     E1 = E1 + 0.5D0 * rho_a0(ii,ii) * FockHF_a0(1,ii,ii)
     E2 = E2 + 0.5D0 * rho_a0(ii,ii) * FockHF_a0(2,ii,ii)
     E3 = E3 + 0.5D0 * rho_a0(ii,ii) * FockHF_a0(3,ii,ii)
     do jj=1,ii-1
       ! Exact HF energy 
       E1 = E1 + 0.5D0 * rho_a0(ii,jj) * FockHF_a0(1,ii,jj)
       E1 = E1 + 0.5D0 * rho_a0(jj,ii) * FockHF_a0(1,jj,ii)
       ! Short HF energy 
       E2 = E2 + 0.5D0 * rho_a0(ii,jj) * FockHF_a0(2,ii,jj)
       E2 = E2 + 0.5D0 * rho_a0(jj,ii) * FockHF_a0(2,jj,ii)
       ! Long HF energy 
       E3 = E3 + 0.5D0 * rho_a0(ii,jj) * FockHF_a0(3,ii,jj)
       E3 = E3 + 0.5D0 * rho_a0(jj,ii) * FockHF_a0(3,jj,ii)
     enddo
   enddo

   if ( OPEN ) then
      do ii=1,M
        E1 = E1 + 0.5D0 * rho_b0(ii,ii) * FockHF_b0(1,ii,ii)
        E2 = E2 + 0.5D0 * rho_b0(ii,ii) * FockHF_b0(2,ii,ii)
        E3 = E3 + 0.5D0 * rho_b0(ii,ii) * FockHF_b0(3,ii,ii)
        do jj=1,ii-1
          ! Exact HF energy 
          E1 = E1 + 0.5D0 * rho_b0(ii,jj) * FockHF_b0(1,ii,jj)
          E1 = E1 + 0.5D0 * rho_b0(jj,ii) * FockHF_b0(1,jj,ii)
          ! Short HF energy 
          E2 = E2 + 0.5D0 * rho_b0(ii,jj) * FockHF_b0(2,ii,jj)
          E2 = E2 + 0.5D0 * rho_b0(jj,ii) * FockHF_b0(2,jj,ii)
          ! Long HF energy 
          E3 = E3 + 0.5D0 * rho_b0(ii,jj) * FockHF_b0(3,ii,jj)
          E3 = E3 + 0.5D0 * rho_b0(jj,ii) * FockHF_b0(3,jj,ii)
        enddo
      enddo
   endif
   E1 = E1 * (-HF_fac(1))
   E2 = E2 * (-HF_fac(2))
   E3 = E3 * (-HF_fac(3))
   call g2g_timer_sum_pause("All Exact Exchange Energy")

end subroutine exact_energies

subroutine exact_exchange_forces(ffx,Pmat_vec,M,natom)
use extern_functional_data, only: HF, HF_fac
   implicit none

   integer, intent(in)  :: M, natom
   LIODBLE, intent(in)  :: Pmat_vec(:)
   LIODBLE, intent(out) :: ffx(natom,3)

   integer :: ii
   logical :: exchange_forces
   LIODBLE, dimension(:,:), allocatable :: rho, frc
    
   ffx = 0.0d0
   exchange_forces = .false.
   do ii=1,3
      if ( HF(ii) == 1 ) exchange_forces = .true.
   enddo
   if (.not. exchange_forces) return

   call g2g_timer_sum_start("All Exact Exchange Gradients")
   allocate(rho(M,M),frc(natom,3))
   call spunpack_rho('L',M,Pmat_vec,rho)

   ! Exact Full HF
   if ( HF(1) == 1 ) then
      frc = 0.0d0
      call g2g_exact_exchange_gradient(rho,frc,1)
      ffx = HF_fac(1) * frc
   endif

   ! Exact Short HF
   if ( HF(2) == 1 ) then
      frc = 0.0d0
      call g2g_exact_exchange_gradient(rho,frc,2)
      ffx = ffx + HF_fac(2) * frc
   endif

   ! Exact Long HF
   if ( HF(3) == 1 ) then
      frc = 0.0d0
      call g2g_exact_exchange_gradient(rho,frc,3)
      ffx = ffx + HF_fac(3) * frc
   endif
   deallocate(rho,frc)
   call g2g_timer_sum_stop("All Exact Exchange Gradients")

end subroutine exact_exchange_forces

subroutine excited_gradients(rhoG,DiffExc,Xmat,fEE,M,natom)
use extern_functional_data, only: HF
   implicit none

   integer, intent(in) :: M, natom
   LIODBLE, intent(in) :: rhoG(M,M), DiffExc(M,M), Xmat(M,M)
   LIODBLE, intent(out) :: fEE(natom,3)

   integer :: ii
   logical :: exact 

   ! Need exact HF gradients ?
   fEE = 0.0d0
   exact = .false.
   do ii=1,3
      if ( HF(ii) == 1 ) exact = .true.
   enddo
   if ( .not. exact ) return

   ! This routine calculates all exact HF gradient terms
   call g2g_exacgrad_excited(rhoG,DiffExc,Xmat,fEE)
end subroutine excited_gradients

end module extern_functional_subs
