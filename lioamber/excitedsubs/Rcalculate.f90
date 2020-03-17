subroutine RCalculate(FXAB,FXIJ,FTIA,GXCIA,X,Rvec,Qvec,NCO,Nvirt,Ndim)
   implicit none

   integer, intent(in) :: NCO, Nvirt, Ndim
   LIODBLE, intent(in) :: FXAB(Nvirt,Nvirt), FXIJ(NCO,NCO)
   LIODBLE, intent(in) :: FTIA(NCO,Nvirt), GXCIA(NCO,Nvirt)
   LIODBLE, intent(in) :: X(Ndim)
   LIODBLE, intent(out) :: Rvec(Ndim), Qvec(Ndim)

   integer :: ii, aa, bb, jj, posf, pos1, NCOc
   LIODBLE :: temp1, temp2

   temp1 = 0.0D0; temp2 = 0.0D0
   NCOc = NCO + 1

   do ii=1,NCO
   do aa=1,Nvirt
      posf = (ii-1) * Nvirt + aa
      ! VIRTUAL PART
      do bb=1,Nvirt
         pos1 = (ii-1) * Nvirt + bb
         temp1 = temp1 + X(pos1) * FXAB(aa,bb)
      enddo
      ! OCC PART
      do jj=1,NCO
         pos1 = (jj-1) * Nvirt + aa
         temp2 = temp2 + X(pos1) * FXIJ(NCOc-ii,NCOc-jj)
      enddo
      Qvec(posf) = temp2
      ! OBTAIN RIA IN VECTOR FORM
      Rvec(posf) = temp2 - (temp1 + FTIA(NCOc-ii,aa) + 2.0D0*GXCIA(NCOc-ii,aa))
      temp1 = 0.0D0; temp2 = 0.0D0
   enddo
   enddo
end subroutine RCalculate
