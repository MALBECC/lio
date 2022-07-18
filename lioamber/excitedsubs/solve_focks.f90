subroutine solve_focks(MatCoef,tvecMO,AX,M,Mlr,NCO,Nvirt,Ndim,&
                       max_subs, vec_dim,Subdim,first_vec)
use excited_data, only: fittExcited
use extern_functional_data, only: need_HF
   implicit none

   integer, intent(in) :: M, Mlr, NCO, Nvirt, Ndim, max_subs, &
                          vec_dim, Subdim, first_vec
   LIODBLE, intent(in) :: MatCoef(M,Mlr), tvecMO(Ndim,max_subs)
   LIODBLE, intent(inout) :: AX(Ndim,max_subs)

   integer :: ivec
   logical :: is_calc
   LIODBLE, allocatable :: tmatMO(:,:), tmatAO(:,:,:)
   LIODBLE, allocatable :: F2e(:,:,:), Fxc(:,:), Ftot(:,:)

   ! Allocate Memory
   allocate(tmatMO(Mlr,Mlr),tmatAO(M,M,vec_dim))
   allocate(F2e(M,M,vec_dim),Fxc(M,M),Ftot(M,M))

   ! Obtain Transition density of all vectors and change to 
   ! AO basis
   call g2g_timer_start("Change Vector(MO) to Mat(AO)")
   do ivec = 1, vec_dim
      call vecMOtomatMO(tvecMO,tmatMO,Mlr,NCO,Nvirt,&
                        Subdim,first_vec,ivec,Ndim)

      ! Change Basis: MO -> AO
      call matMOtomatAO(tmatMO,tmatAO(:,:,ivec),MatCoef,M,Mlr,.true.)
   enddo
   call g2g_timer_stop("Change Vector(MO) to Mat(AO)")

   ! Calculate 2E Integrals of all inner vectors
   is_calc = .false.
   if ( .not. fittExcited ) then
      call g2g_timer_start("Fock 2e LR")
      call g2g_calculate2e(TmatAO,F2e,vec_dim)
      call g2g_timer_stop("Fock 2e LR")
      is_calc = .true.
   endif 

   ! Start inner vectors loop
   do ivec = 1, vec_dim
      ! Calculate 2E Integrals
      if ( .not. is_calc ) then
         if ( fittExcited .and. (.not. need_HF) ) then
            call g2g_timer_start("Fock 2e LR")
            call calc2eFITT(tmatAO(:,:,ivec),F2e(:,:,ivec),M)
            call g2g_timer_stop("Fock 2e LR")
         else
            print*, "Error in 2 Electron Repulsion Integrals"
            print*, "Check HF in the functional and fittExcited"
            stop
         endif
      endif

      ! Calculate XC Integrals
      call g2g_timer_start("Fock XC LR")
      Fxc = 0.0d0
      call g2g_calculateXC(tmatAO(:,:,ivec),Fxc,2)
      call g2g_timer_stop("Fock XC LR")
 
      ! Total Fock
      Ftot = F2e(:,:,ivec) + Fxc

      ! AX = CT * F * C
      call MtoIANV(Ftot,MatCoef,AX,M,Mlr,NCO,Ndim,Subdim,first_vec,ivec)
   enddo

   deallocate(tmatMO,tmatAO,F2e,Fxc,Ftot)
end subroutine solve_focks
