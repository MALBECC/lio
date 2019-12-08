subroutine solve_focks(MatCoef,tvecMO,AX,M,NCO,Nvirt,Ndim,&
                       max_subs,nstates,vec_dim,Subdim,first_vec)
use excited_data, only: fittExcited
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
      if ( .not. fittExcited ) then
         call g2g_calculate2e(tmatAO,F2e)
      else
         print*, "FITT NOT IMPLEMENTED YED"
         stop
      endif

      ! Calculate XC Integrals
      Fxc = 0.0d0
      call g2g_calculateXC(tmatAO,Fxc)
 
      ! Total Fock
      Ftot = F2e + Fxc

      ! AX = CT * F * C
      call MtoIANV(Ftot,MatCoef,AX,M,NCO,Ndim,Subdim,first_vec,ivec)
   enddo

   deallocate(tmatMO,tmatAO,F2e,Fxc,Ftot)
end subroutine solve_focks
