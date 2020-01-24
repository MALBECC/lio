subroutine matMOtomatAO(MatMO,MatAO,C,M,Mlr,transp)
use excited_data, only: Coef_trans

   implicit none

   integer, intent(in) :: M, Mlr
   logical, intent(in) :: transp
   double precision, intent(in) :: MatMO(Mlr,Mlr), C(M,Mlr)
   double precision, intent(out) :: MatAO(M,M)

   double precision, dimension(:,:), allocatable :: scratch

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
