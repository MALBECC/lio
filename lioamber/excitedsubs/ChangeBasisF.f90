subroutine ChangeBasisF(FX,FT,Gxc,FXAB,FXIJ,FTIA,GXCIA,M,Nvirt,NCO)
use excited_data          , only: Cocc, Cocc_trans, Cvir, Cvir_trans
use extern_functional_data, only: HF

   implicit none

   integer, intent(in) :: M, Nvirt, NCO
   LIODBLE, intent(in)  :: Gxc(M,M), FX(M,M), FT(M,M)
   LIODBLE, intent(out) :: FXAB(Nvirt,Nvirt), FXIJ(NCO,NCO)
   LIODBLE, intent(out) :: FTIA(NCO,Nvirt), GXCIA(NCO,Nvirt)

   character(len=1) :: code
   integer :: i
   LIODBLE, dimension(:,:), allocatable :: scratch

   code = 'N'
   do i=1,3
      if ( HF(i)==1 ) code = 'T'
   enddo

!  FORM FX IN BASIS VIRT X VIRT
   allocate(scratch(M,Nvirt))
   call dgemm(code,'N',M,Nvirt,M,1.0d0,FX,M,Cvir,M,0.0d0,scratch,M)
   call dgemm('N','N',Nvirt,Nvirt,M,1.0d0,Cvir_trans,Nvirt,scratch,&
              M,0.0d0,FXAB,Nvirt)
   deallocate(scratch)

!  FORM FX IN BASIS OCC X OCC
   allocate(scratch(M,NCO))
   call dgemm(code,'N',M,NCO,M,1.0d0,FX,M,Cocc,M,0.0d0,scratch,M)
   call dgemm('N','N',NCO,NCO,M,1.0d0,Cocc_trans,NCO,scratch,M,&
              0.0d0,FXIJ,NCO)
   deallocate(scratch)

!  FORM FT IN BASIS OCC X VIR
   allocate(scratch(M,Nvirt))
   call dgemm(code,'N',M,Nvirt,M,1.0d0,FT,M,Cvir,M,0.0d0,scratch,M)
   call dgemm('N','N',NCO,Nvirt,M,1.0d0,Cocc_trans,NCO,scratch,M,&
              0.0d0,FTIA,NCO)
   deallocate(scratch)

!  FORM GXC IN BASIS OCC X VIR
   allocate(scratch(M,Nvirt))
   call dgemm('N','N',M,Nvirt,M,1.0d0,Gxc,M,Cvir,M,0.0d0,scratch,M)   
   call dgemm('N','N',NCO,Nvirt,M,1.0d0,Cocc_trans,NCO,scratch,M, &
             0.0d0,GXCIA,NCO)
   deallocate(scratch)

end subroutine ChangeBasisF
