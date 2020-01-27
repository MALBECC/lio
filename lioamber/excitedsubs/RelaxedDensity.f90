subroutine RelaxedDensity(Xexc,Eexc,C_scf,E_scf,M,Mlr,Nvirt,NCO,Ndim,Nstat)
! This routine generates the relaxed density matrix of the root 
! Excited State
use excited_data, only: root
   implicit none

   integer, intent(in) :: M, Mlr, Nvirt, NCO, Ndim, Nstat
   double precision, intent(in) :: Xexc(Ndim,Nstat), Eexc(Nstat)
   double precision, intent(in) :: C_scf(M,Mlr), E_scf(Mlr)

   double precision, allocatable :: Xlr(:), Punr(:,:), Trans(:,:)
   double precision, allocatable :: Zvec(:), Qvec(:), Gxc(:,:)
   double precision, allocatable :: rhoEXC(:,:), Pdif(:,:)

   if ( root == 0 ) return
   if ( root > Nstat ) then
      print*, "The excited state that you want wasn't &
              & calculated in Linear Response"
      print*, "Please, setup root <= nstates"
      stop
   endif

   write(*,"(1X,A,1X,I2)") "FORM RELAXED DENSITY MATRIX FOR EXCITED STATE:", root
   write(*,*) ""

   allocate(Xlr(Ndim)); Xlr = Xexc(:,root) / dsqrt(2.0d0)
   allocate(Punr(M,M),Trans(M,M)); Punr = 0.0d0; Trans = 0.0d0
   ! Xlr = Transition Density in vector form, in MO basis
   ! Punr = Unrelaxed Difference Density Matrix, in AO basis
   ! Trans = Transition Density Matrix, in AO basis
   call GenerateDensities(Xlr,C_scf,Punr,Trans,M,Mlr,Ndim,NCO,Nvirt)
   allocate(Zvec(Ndim),Qvec(Ndim),Gxc(M,M))
   Zvec = 0.0d0; Qvec = 0.0d0; Gxc = 0.0d0
   call Zvector(C_scf,E_scf,Xlr,Punr,Trans,Zvec,Qvec,Gxc,NCO,&
                M, Mlr, Ndim, Nvirt)

   ! Obtain Densities
   allocate(rhoEXC(M,M),Pdif(M,M)); rhoEXC = 0.0d0; Pdif = 0.0d0
   call ObtainDens(Zvec,Punr,C_scf,rhoEXC,Pdif,M,Mlr,NCO,Ndim,Nvirt)

end subroutine

subroutine ObtainDens(Z,Rho_urel,C,Rho_exc,Rel_diff,M,Mlr,NCO,N,Nvirt)
use garcha_mod, only: Pmat_vec
   implicit none

   integer, intent(in) :: M, Mlr, NCO, N, Nvirt
   double precision, intent(in) :: Z(N), C(M,Mlr), Rho_urel(M,M)
   double precision, intent(out) :: Rho_exc(M,M), Rel_diff(M,M)

   integer :: i, j, NCOc, pos
   double precision, dimension(:,:), allocatable :: Zmo, Zao, Rho_fund

   ! Extract rho GS
   allocate(Rho_fund(M,M))
   call spunpack_rho('L', M, Pmat_vec, Rho_fund)

   ! Convert Z in AO Basis
   allocate(Zmo(M,M)); Zmo = 0.0D0
   NCOc = NCO + 1
   do i=1,NCO
   do j=1,Nvirt
      pos = (i - 1) * Nvirt + j
      Zmo(NCOc-i,NCO+j) = Z(pos)
   enddo
   enddo
   allocate(Zao(M,M))
   call matMOtomatAO(Zmo,Zao,C,M,Mlr,.false.)
   deallocate(Zmo)

   Rel_diff = Rho_urel + Zao
   deallocate(Zao)

   Rho_exc = Rho_fund + Rel_diff + transpose(Rel_diff)
   
   deallocate(Rho_fund)
end subroutine ObtainDens
