subroutine addInt(A,Ene,Vec,Ndim,Mlr,Sdim,NCO,Nvec,V1)
   implicit none

   integer, intent(in) :: Ndim, Mlr, NCO, Nvec, V1, Sdim
   double precision, intent(in) :: Ene(Mlr), Vec(Ndim,Sdim)
   double precision, intent(inout) :: A(Ndim,Sdim)

   integer :: iv, i, j, Nvirt, NCOc, row

   Nvirt = Mlr - NCO
   NCOc = NCO + 1
   do iv=1,Nvec
     do i=1,NCO
     do j=NCOc,Mlr
       row = (i-1) * Nvirt + (j-NCO)
       A(row,V1+iv) = A(row,V1+iv) + (Ene(j) - Ene(NCOc-i)) * Vec(row,V1+iv)
     enddo
     enddo
   enddo
end subroutine addInt

