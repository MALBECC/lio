!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module mask_ecp
   implicit none
   contains


!------------------------------------------------------------------------------!
   subroutine ECP_init()

      use ECP_mod, only : ecpmode, FOCK_ECP_read, FOCK_ECP_write
      implicit none

      if (.not.ecpmode) return

      if (FOCK_ECP_read) then
!        alocatea variables comunes y las lee del archivo ECP_restart
         call intECP(0)
      else
!        intECP(1) alocatea variables, calcula variables comunes, y calcula
!        terminos de 1 centro, mientras que intECP(2/3) calcula los términos
!        de 2 y 3 centros respectivamente
         call g2g_timer_start('ECP Routines')
         call intECP(1)
         call intECP(2)
         call intECP(3)
         call g2g_timer_stop('ECP Routines')
      end if

      if (FOCK_ECP_write) then
         call WRITE_ECP()
      end if

      call WRITE_POST(1)
   end subroutine ECP_init


!------------------------------------------------------------------------------!
   subroutine ECP_fock( Nvec, fockvec )

      use ECP_mod, only : ecpmode, term1e, VAAA, VAAB, VBAC
      implicit none
      integer, intent(in)    :: Nvec 
      real*8 , intent(inout) :: fockvec(Nvec)
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
      real*8 , intent(in)    :: densvec(Nvec)
      real*8 , intent(inout) :: E_ecp
      real*8 , intent(inout) :: E_mod
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
