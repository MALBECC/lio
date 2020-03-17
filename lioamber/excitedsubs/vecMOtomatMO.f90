subroutine vecMOtomatMO(vec,mat,Mlr,NCO,Nvirt,Sdim,V1,iv,dim)
   implicit none

   integer, intent(in) :: Mlr, NCO, Nvirt, Sdim, V1, iv, dim
   LIODBLE, intent(in)  :: vec(dim,Sdim)
   LIODBLE, intent(out) :: mat(Mlr,Mlr)

   integer :: row, col, NCOc

   NCOc = NCO - 1
   mat = 0.0D0

   do row=0,NCOc
     do col=1,Nvirt
        mat(NCOc-row+1,NCO+col) = vec(row*Nvirt+col,V1+iv)
     enddo
   enddo
end subroutine vecMOtomatMO
