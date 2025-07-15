subroutine solve_focks(MatCoef,tvecMO,AX,M,Mlr,NCO,Nvirt,Ndim,&
                       max_subs, vec_dim,Subdim,first_vec)
use excited_data, only: fittExcited, Coef_trans, Cocc_trans
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
      call matMOtomatAO(tmatMO,tmatAO(:,:,ivec),MatCoef,Coef_trans,M,Mlr,.true.)
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
      call MtoIANV(Ftot,MatCoef,Cocc_trans,AX,M,Mlr,NCO,Ndim,Subdim,first_vec,ivec)
   enddo

   deallocate(tmatMO,tmatAO,F2e,Fxc,Ftot)
end subroutine solve_focks

! OPEN SHELL SOLVE FOCKS
subroutine open_solve_focks(MatCoefA,tvecMOA,MatCoefB,tvecMOB,AX_a,AX_b, &
                            M,Mlr,NCOA,NvirtA,NCOB,NvirtB,Ndim, &
                            max_subs,vec_dim,Subdim,first_vec)
use excited_data, only: fittExcited, Coef_trans, Coef_transB, Cocc_trans, Cocc_transB
use extern_functional_data, only: need_HF
   implicit none

   integer, intent(in) :: M, Mlr, NCOA, NvirtA, NCOB, NvirtB, Ndim, max_subs, &
                          vec_dim, Subdim, first_vec
   LIODBLE, intent(in) :: MatCoefA(M,Mlr), tvecMOA(Ndim,max_subs)
   LIODBLE, intent(in) :: MatCoefB(M,Mlr), tvecMOB(Ndim,max_subs)
   LIODBLE, intent(inout) :: AX_a(Ndim,max_subs), AX_b(Ndim,max_subs)

   integer :: ivec
   logical :: is_calc
   LIODBLE, allocatable :: tmatMOa(:,:), tmatAOa(:,:,:), tmatMOb(:,:), tmatAOb(:,:,:)
   LIODBLE, allocatable :: F2eA(:,:,:), FxcA(:,:), FtotA(:,:)
   LIODBLE, allocatable :: F2eB(:,:,:), FxcB(:,:), FtotB(:,:)

   ! Allocate Memory
   allocate(tmatMOa(Mlr,Mlr),tmatAOa(M,M,vec_dim),tmatMOb(Mlr,Mlr),tmatAOb(M,M,vec_dim))
   allocate(F2eA(M,M,vec_dim),FxcA(M,M),FtotA(M,M),F2eB(M,M,vec_dim),FxcB(M,M),FtotB(M,M))

   ! Obtain Transition density of all vectors and change to 
   ! AO basis for both spin
   call g2g_timer_start("Change Vector(MO) to Mat(AO)")
   do ivec = 1, vec_dim 
      ! For alpha
      call vecMOtomatMO(tvecMOA,tmatMOa,Mlr,NCOA,NvirtA,&
                        Subdim,first_vec,ivec,Ndim)
      call matMOtomatAO(tmatMOa,tmatAOa(:,:,ivec),MatCoefA,Coef_trans,M,Mlr,.true.)

      ! For beta
      call vecMOtomatMO(tvecMOB,tmatMOb,Mlr,NCOB,NvirtB,&
                        Subdim,first_vec,ivec,Ndim)
      call matMOtomatAO(tmatMOb,tmatAOb(:,:,ivec),MatCoefB,Coef_transB,M,Mlr,.true.)
   enddo
   call g2g_timer_stop("Change Vector(MO) to Mat(AO)")

   ! Calculate 2E Integrals of all inner vectors
   call g2g_timer_start("Fock 2e LR")
   call g2g_calculate2e_open(tmatAOa,tmatAOb,F2eA,F2eB,vec_dim)
   call g2g_timer_stop("Fock 2e LR")

   ! Calculate XC Integrals
   do ivec = 1, vec_dim
      call g2g_timer_start("Fock XC LR")
      FxcA = 0.0d0; FxcB = 0.0d0
      call g2g_open_calculateXC(tmatAOa(:,:,ivec),tmatAOb(:,:,ivec),FxcA,FxcB,2)
      call g2g_timer_stop("Fock XC LR")

      ! Total Fock
      FtotA = F2eA(:,:,ivec) + FxcA
      FtotB = F2eB(:,:,ivec) + FxcB

      ! AX_a,b = CT * F * C
      call MtoIANV(FtotA,MatCoefA,Cocc_trans,AX_a,M,Mlr,NCOA,Ndim,Subdim,first_vec,ivec)
      call MtoIANV(FtotB,MatCoefB,Cocc_transB,AX_b,M,Mlr,NCOB,Ndim,Subdim,first_vec,ivec)
   enddo
   deallocate(tmatMOa,tmatAOa,F2eA,FxcA,FtotA)
   deallocate(tmatMOb,tmatAOb,F2eB,FxcB,FtotB)
end subroutine open_solve_focks

