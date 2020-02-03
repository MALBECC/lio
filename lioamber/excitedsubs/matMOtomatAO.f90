subroutine matMOtomatAO(MatMO,MatAO,C,M,transp)
use excited_data, only: Coef_trans

   implicit none

   integer, intent(in) :: M
   logical, intent(in) :: transp
   double precision, intent(in) :: MatMO(M,M), C(M,M)
   double precision, intent(out) :: MatAO(M,M)

   double precision, dimension(:,:), allocatable :: scratch

!  CHANGE BASIS MO -> AO
   allocate(scratch(M,M)); MatAO = 0.0d0
   call dgemm('N','N',M,M,M,1.0d0,MatMO,M,Coef_trans,M,0.0d0,scratch,M)
   call dgemm('N','N',M,M,M,1.0d0,C,M,scratch,M,0.0d0,MatAO,M)

!  WE TRANSPOSE MATRIX FOR USE IT IN C 
!  TODO: DELETE THIS
   if ( transp ) then
      scratch = MatAO(:,:)
      MatAO   = transpose(scratch)
   endif

   deallocate(scratch)
end subroutine matMOtomatAO
