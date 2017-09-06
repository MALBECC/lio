!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module dftb_subs
   implicit none

contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dftb_init()
end subroutine dftb_init
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dftb_final()
end subroutine dftb_final
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine find_neighbors(M,Nuc,natom, nend, Iend)

   implicit none
   
   integer              ::  i, j
   integer, intent(in)  ::  M, natom
   integer, intent(in)  ::  Nuc(M)
   integer, intent(in) ::  nend
   integer, intent(out) :: Iend(2 , nend)

!   nend=0

!   do i=1, M
!      if (Nuc(i)==1.or.Nuc(i)==2.or.Nuc(i)==natom-1.or.Nuc(i)==natom) then
!         nend = 1 + nend
!      endif
!   end do

!   allocate (iend(2, nend)) 
   
   j=0
   do i=1, M
      if (Nuc(i)==1.or.Nuc(i)==2.or.Nuc(i)==natom-1.or.Nuc(i)==natom) then
         j=j+1
         Iend(1,j) = Nuc(i)
         Iend(2,j) = i
      end if
   end do   

end subroutine find_neighbors

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine build_chimera (MTB, M, MTBQM, alfaTB, betaTB, gammaTB, fock,      &
                          chimerafock, Iend, nend, natom)

   implicit none

   integer, intent(in)  :: MTB, M, MTBQM, nend, natom !MTB debe tener la dimension de
!los atomos de los bias, sino contar el que esta dentro de fock
   integer, intent(in)  :: Iend(2, nend)
   real*8, intent(in)   :: alfaTB, betaTB, gammaTB
   real*8, intent(inout)   :: fock (M, M)
   integer              :: i, j

!#ifdef TD_SIMPLE
   real*8, intent(out)  :: chimerafock (MTBQM, MTBQM) !temporal dimensions
!#else
!   real*16, intent(out) :: chimerafock (2*MTB+M, 2*MTB+M)
!#endif

   chimerafock(:,:) = 0.0D0

   do i=1, nend
      if (Iend(1, i) == 1) then

         fock(Iend(2, i),:) = 0.0D0
         fock(:,Iend(2, i)) = 0.0D0
         fock(Iend(2, i),Iend(2, i))= alfaTB
         
         chimerafock(Iend(2,i)+MTB,MTB) = betaTB
         chimerafock(MTB, Iend(2,i)+MTB) = betaTB
 
         do j=1, nend
           if (Iend(1, j) == 2) then
              fock(Iend(2, i), Iend(2, j))=gammaTB
              fock(Iend(2, j), Iend(2, i))=gammaTB
           end if
         end do

      else if (Iend(1, i) == natom) then

         fock(Iend(2, i),:) = 0.0D0
         fock(:,Iend(2, i)) = 0.0D0
         fock(Iend(2, i),Iend(2, i))= alfaTB

         chimerafock(Iend(2,i)+MTB,MTB+M+1) = betaTB
         chimerafock(MTB+M+1, Iend(2,i)+MTB) = betaTB

         do j=1, nend
           if (Iend(1, j) == natom-1) then
              fock(Iend(2, i), Iend(2, j))=gammaTB
              fock(Iend(2, j), Iend(2, i))=gammaTB
           end if
         end do
      
      end if
   end do


   do i=1,MTB
      chimerafock(i,i) = alfaTB
      chimerafock(MTB+M+i, MTB+M+i)= alfaTB
      
      if (i<MTB) then

         chimerafock(i,i+1) = betaTB
         chimerafock(i+1,i) = betaTB 
         chimerafock(MTB+M+i, MTB+M+i+1)= betaTB
         chimerafock(MTB+M+i+1, MTB+M+i)= betaTB

      end if
   end do


   chimerafock(MTB+1:MTB+M, MTB+1:MTB+M) = fock(:,:)

end subroutine build_chimera

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine disarm_chimera (MTB, M, fock, chimerafock)

   implicit none
   integer, intent(in)  :: MTB, M
   real*8, intent(out) :: fock (M, M)
   integer            :: i, j

!#ifdef TD_SIMPLE
   real*8, intent(in)  :: chimerafock (2*MTB+M, 2*MTB+M) !temporal dimensions
!#else
!   real*16, intent(in) :: chimerafock (2*MTB+M, 2*MTB+M)
!#endif

   fock(:,:) = chimerafock(MTB+1:MTB+M, MTB+1:MTB+M)

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module dftb_subs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
