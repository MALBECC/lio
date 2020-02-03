subroutine addInt(A,Ene,Vec,Ndim,M,Sdim,NCO,Nvec,V1)
   implicit none

   integer, intent(in) :: Ndim, M, NCO, Nvec, V1, Sdim
   double precision, intent(in) :: Ene(M), Vec(Ndim,Sdim)
   double precision, intent(inout) :: A(Ndim,Sdim)

   integer :: iv, i, j, Nvirt, NCOc, row

   Nvirt = M - NCO
   NCOc = NCO + 1
   do iv=1,Nvec
     do i=1,NCO
     do j=NCOc,M
       row = (i-1) * Nvirt + (j-NCO)
       A(row,V1+iv) = A(row,V1+iv) + (Ene(j) - Ene(NCOc-i)) * Vec(row,V1+iv)
     enddo
     enddo
   enddo
end subroutine addInt

