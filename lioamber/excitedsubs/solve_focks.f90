subroutine solve_focks(MatCoef,tvecMO,AX,M,NCO,Nvirt,Ndim,&
                       max_subs,nstates,vec_dim,Subdim,first_vec)
use excited_data, only: fittExcited
use garcha_mod,   only: PBE0
   implicit none

   integer, intent(in) :: M, NCO, Nvirt, Ndim, max_subs, nstates, vec_dim, &
                          Subdim, first_vec
   double precision, intent(in) :: MatCoef(M,M), tvecMO(Ndim,max_subs)
   double precision, intent(inout) :: AX(Ndim,max_subs)

   integer :: ivec
   double precision, allocatable :: tmatMO(:,:), tmatAO(:,:)
   double precision, allocatable :: F2e(:,:), Fxc(:,:), Ftot(:,:)

   ! Allocate Memory
   allocate(tmatMO(M,M),tmatAO(M,M))
   allocate(F2e(M,M),Fxc(M,M),Ftot(M,M))

   ! Start loop of inner vectors
   do ivec = 1, vec_dim
      ! Change Vectors to Matrix in OM
      call vecMOtomatMO(tvecMO,tmatMO,M,NCO,Nvirt,&
                        Subdim,first_vec,ivec,Ndim)

      ! Change Basis: MO -> AO
      call matMOtomatAO(tmatMO,tmatAO,MatCoef,M,.true.)
      
      ! Calculate 2E Integrals
      call g2g_timer_start("Fock 2e LR")
      if ( .not. fittExcited ) then
         call g2g_calculate2e(tmatAO,F2e)
      elseif ( fittExcited .and. (.not. PBE0) ) then
         call calc2eFITT(tmatAO,F2e,M)
      else
         print*, "Error in 2 Electron Repulsion Integrals"
         print*, "Check PBE0 and fittExcited"
         stop
      endif
      call g2g_timer_stop("Fock 2e LR")

      ! Calculate XC Integrals
      call g2g_timer_start("Fock XC LR")
      Fxc = 0.0d0
      call g2g_calculateXC(tmatAO,Fxc)
      call g2g_timer_stop("Fock XC LR")
 
      ! Total Fock
      Ftot = F2e + Fxc

      ! AX = CT * F * C
      call MtoIANV(Ftot,MatCoef,AX,M,NCO,Ndim,Subdim,first_vec,ivec)
   enddo

   deallocate(tmatMO,tmatAO,F2e,Fxc,Ftot)
end subroutine solve_focks
