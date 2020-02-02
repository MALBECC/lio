subroutine HVgradcalc(rho,for,M,natom,nuclear)
use garcha_mod, only: d, r, Iz, ntatom
use faint_cpu , only: int1G
   implicit none

   integer, intent(in) :: M, natom
   logical, intent(in) :: nuclear
   double precision, intent(in)  :: rho(M,M)
   double precision, intent(out) :: for(natom,3)

   integer :: ii, jj, ind, MM
   double precision, allocatable :: rho_vec(:)

   MM = M * (M + 1) / 2
   allocate(rho_vec(MM))
   ind = 1
   do ii=1,M
      rho_vec(ind) = rho(ii,ii)
      ind = ind + 1
      do jj=ii+1,M
         rho_vec(ind) = rho(ii,jj) * 2.0d0
         ind = ind + 1
      enddo
   enddo

   ! GRADIENTS
   for = 0.0d0
   call int1G(for,rho_vec,d,r,Iz,natom,ntatom,nuclear,.true.)
   deallocate(rho_vec)
end subroutine HVgradcalc
