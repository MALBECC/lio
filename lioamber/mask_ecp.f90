#include "datatypes/datatypes.fh"
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module mask_ecp
   implicit none
   contains

!------------------------------------------------------------------------------!
   subroutine ECP_init()

      use ECP_mod, only : ecpmode, FOCK_ECP_read, FOCK_ECP_write, &
#ifdef FULL_CHECKS
      inf_Q, NAN_Q, inf_Q2, NAN_Q2, &
#endif
      FOCK_ECP_write, first_steep

      use faint_cpu  , only: intECPG

      implicit none

      if (.not.ecpmode) return

      if (FOCK_ECP_read) then
         call generalECP(0) ! alocatea variables comunes y las lee del archivo ECP_restart, solo para calculos sin gradiente
      else
         call g2g_timer_start('ECP Routines')
         if (first_steep) then
            call generalECP(1) ! alocatea variables, calcula variables comunes, y calcula terminos de 1 centro.
            first_steep=.false.
         endif
         call generalECP(5) ! calcula terminos de 2 y 3 centros y sus derivadas.
         call g2g_timer_stop('ECP Routines')
      end if

      if (FOCK_ECP_write) then
         call WRITE_ECP()
      end if

      call WRITE_POST(1)

#ifdef FULL_CHECKS
      write(*,*) "Reporting error in radial integrals"
      write(*,*) "infty in Q: ", inf_Q ,"NaN in Q: ", NAN_Q
      write(*,*) "infty in Q2: ", inf_Q2 , "NaN in Q2: ", NAN_Q2
#endif

   end subroutine ECP_init


!------------------------------------------------------------------------------!
   subroutine ECP_fock( Nvec, fockvec )

      use ECP_mod, only : ecpmode, term1e, VAAA, VAAB, VBAC
      implicit none
      integer, intent(in)    :: Nvec 
      LIODBLE , intent(inout) :: fockvec(Nvec)
      integer                :: kk

      if (.not.ecpmode) return

!     backup of 1e terms and adds ECP to fock
      write(*,*) "Modifying Fock Matrix with ECP terms"
      do kk = 1, Nvec
         term1e(kk)  = fockvec(kk)
         fockvec(kk) = fockvec(kk) + VAAA(kk) + VAAB(kk) + VBAC(kk)
      end do

   end subroutine ECP_fock


!------------------------------------------------------------------------------!
   subroutine ECP_energy( Nvec, densvec, E_ecp, E_mod )

      use ECP_mod, only : ecpmode, VAAA, VAAB, VBAC
      implicit none
      integer, intent(in)    :: Nvec
      LIODBLE , intent(in)    :: densvec(Nvec)
      LIODBLE , intent(inout) :: E_ecp
      LIODBLE , intent(inout) :: E_mod
      integer                :: kk

      if (.not.ecpmode) return

      E_ecp = 0.d0
      do kk =1, Nvec
         E_ecp = E_ecp + densvec(kk) * (VAAA(kk)+VAAB(kk)+VBAC(kk))
      end do
      E_mod = E_mod - E_ecp

   end subroutine ECP_energy


end module mask_ecp
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
