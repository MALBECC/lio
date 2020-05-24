subroutine getcubegen_excited(Xexc,Trans,Pdif,rhoExc,Coef,M,Ndim,nstates)
use garcha_mod,   only: Pmat_vec
use basis_data,   only: MM
use excited_data, only: cgPe, cgPd, cgPt, cgMO1, cgPg, root
use cubegen_data, only: cube_dens, cube_dens_file, cube_orb, cube_sel, &
                        cube_orb_file
use cubegen,      only: cubegen_write
use SCF_aux,      only: messup_densmat
! Inputs
! Xexc = Eigenvectors of LR
! Trans= Transiton density matrix
! Pdif = Relaxed Difference Density Matrix
! rhoExc= Relaxed Excited State Density Matrix
! Coef = Molecular Orbitals Coefficients
! M    = number of basis functions
! Ndim = number of single excited configurations
! nstates= number of calculated excited states 

! Cubgen
! cgPe = cubegen of rhoExc
! cgPd = cubegen of Pdif
! cgMO1= cubegen of principal MOs involved in each states
! cgPg = cubeg of Ground state density matrix
! cgPt = cubeg of Transition density matrix

! TODO:
! Get the principal contributions to nstates and generates cubegen of these MOs
   implicit none

   integer, intent(in) :: M, Ndim, nstates
   LIODBLE, intent(in) :: Xexc(Ndim,nstates), Trans(M,M), Pdif(M,M)
   LIODBLE, intent(in) :: rhoExc(M,M), Coef(M,M)

   logical :: need_cubegen
   LIODBLE, allocatable :: rho_vec(:), rho_mat(:,:)

   need_cubegen = ( cgPe .or. cgPd .or. cgPt .or. cgMO1 .or. cgPg )
   if (.not. need_cubegen) return

   need_cubegen = ( cgPe .or. cgPd .or. cgPt )
   if ( need_cubegen .and. root == 0 ) then
      print*, "All densities in excited states are computed if root /= 0"
      print*, "Please set this value."
      stop
   endif

   allocate(rho_vec(MM)) ; rho_vec = Pmat_vec
   allocate(rho_mat(M,M)); rho_mat = 0.0d0

   print*, " "
   print*, "Performing cubegen analysis"

   ! Ground Density
   if ( cgPg ) then
      print*, "Cubegen file of Ground State Density Matrix"
      cube_dens = .true.
      cube_dens_file = "GSdensity.cube"

      call cubegen_write(Coef)
   endif

   ! Excited Density
   if ( cgPe ) then
      print*, "Cubegen file of Relaxed Excited State Density Matrix"
      cube_dens = .true.
      cube_dens_file = "ESdensity.cube"

      rho_mat = rhoExc
      call messup_densmat(rho_mat)
      call sprepack('L',M,Pmat_vec,rho_mat)
      call cubegen_write(Coef)
   endif

   ! Difference Matrix
   if ( cgPd ) then
      print*, "Cubegen file of Relaxed Difference Density Matrix"
      cube_dens = .true.
      cube_dens_file = "DIFdensity.cube"

      rho_mat = Pdif + transpose(Pdif)
      call messup_densmat(rho_mat)
      call sprepack('L',M,Pmat_vec,rho_mat)
      call cubegen_write(Coef)
   endif

   ! Transition Matrix
   if ( cgPt ) then
      print*, "Cubegen file of Transition Density Matrix"
      cube_dens = .true.
      cube_dens_file = "TRANdensity.cube"

      rho_mat = Trans + transpose(Trans)
      call messup_densmat(rho_mat)
      call sprepack('L',M,Pmat_vec,rho_mat)
      call cubegen_write(Coef)
   endif

   Pmat_vec = rho_vec
   deallocate(rho_vec,rho_mat)
end subroutine getcubegen_excited
