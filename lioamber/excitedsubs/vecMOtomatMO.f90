subroutine vecMOtomatMO(vec,mat,M,NCO,Nvirt,Sdim,V1,iv,dim)
   implicit none

   integer, intent(in) :: M, NCO, Nvirt, Sdim, V1, iv, dim
   double precision, intent(in)  :: vec(dim,Sdim)
   double precision, intent(out) :: mat(M,M)

   integer :: row, col, NCOc

   NCOc = NCO - 1
   mat = 0.0D0

   do row=0,NCOc
     do col=1,Nvirt
        mat(NCOc-row+1,NCO+col) = vec(row*Nvirt+col,V1+iv)
     enddo
   enddo
end subroutine vecMOtomatMO
