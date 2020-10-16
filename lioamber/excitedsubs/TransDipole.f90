subroutine TransDipole(Tdens,Tdip,M)
   use garcha_mod, only: NCO, Nunp, r, d, Iz, pc
   use properties, only: dipole

   implicit none
   integer, intent(in) :: M
   LIODBLE, intent(inout) :: Tdens(M,M)
   LIODBLE, intent(inout) :: Tdip(3)

   integer :: ii, jj
   LIODBLE, dimension(:), allocatable :: P_dens

   do ii=1,M
   do jj=1,ii-1
      Tdens(ii,jj) = Tdens(ii,jj) + Tdens(jj,ii)
   enddo
   enddo

   allocate(P_dens(M*(M+1)/2)); P_dens = 0.0d0
   call sprepack('L',M,P_dens,Tdens)
   call dipole(Tdip,P_dens, 2*NCO+Nunp, r, d, Iz, pc)
   Tdip = Tdip * 2.0d0 / dsqrt(2.0d0)
   deallocate(P_dens)

end subroutine TransDipole

   


