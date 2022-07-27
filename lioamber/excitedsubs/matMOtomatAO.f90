subroutine matMOtomatAO(MatMO,MatAO,C,Coef_trans,M,Mlr,transp)
   implicit none

   integer, intent(in) :: M, Mlr
   logical, intent(in) :: transp
   LIODBLE, intent(in) :: MatMO(Mlr,Mlr), C(M,Mlr), Coef_trans(Mlr,M)
   LIODBLE, intent(out) :: MatAO(M,M)

   LIODBLE, dimension(:,:), allocatable :: scratch

!  CHANGE BASIS MO -> AO
   allocate(scratch(Mlr,M)); MatAO = 0.0d0
   call dgemm('N','N',Mlr,M,Mlr,1.0d0,MatMO,Mlr,Coef_trans,Mlr,0.0d0,scratch,Mlr)
   call dgemm('N','N',M,M,Mlr,1.0d0,C,M,scratch,Mlr,0.0d0,MatAO,M)

!  WE TRANSPOSE MATRIX FOR USE IT IN C 
!  TODO: DELETE THIS
   if ( transp ) then
      scratch = MatAO(:,:)
      MatAO   = transpose(scratch)
   endif

   deallocate(scratch)
end subroutine matMOtomatAO
