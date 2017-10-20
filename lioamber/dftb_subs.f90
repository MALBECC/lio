!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module dftb_subs
   implicit none

contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dftb_init(M)

   use dftb_data, only: MTB, MDFTB, end_basis, Iend

   implicit none

   integer, intent(in) :: M

   MDFTB=2*MTB+M
   allocate(Iend(2,2*end_basis))
   
end subroutine dftb_init
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine getXY_DFTB(M_in,x_in,y_in,xmat,ymat)

   use dftb_data, only: MTB, MDFTB

   implicit none
   integer, intent(in)  :: M_in
   real*8 , intent(in)  :: x_in(M_in,M_in)
   real*8 , intent(in)  :: y_in(M_in,M_in)
   real*8 , intent(out) :: xmat(MDFTB,MDFTB)
   real*8 , intent(out) :: ymat(MDFTB,MDFTB)
   
   integer              :: ii,jj


   xmat=0.0d0
   ymat=0.0d0

   do ii=1, MTB
      xmat(ii,ii)=1.0d0
      xmat(MTB+M_in+ii,MTB+M_in+ii)=1.0d0

      ymat(ii,ii)=1.0d0
      ymat(MTB+M_in+ii,MTB+M_in+ii)=1.0d0
   end do   

   do jj=1, M_in
   do ii=1, M_in
      xmat(MTB+ii, MTB+jj)=x_in(ii,jj)
      ymat(MTB+ii, MTB+jj)=y_in(ii,jj)
   end do
   end do

end subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine find_neighbors(M_in,Nuc,natom)

   use dftb_data, only:Iend
  
   implicit none
   
   integer              ::  ii, jj
   integer, intent(in)  ::  M_in, natom
   integer, intent(in)  ::  Nuc(M_in)

!   end_basis=0

!   do i=1, M
!      if (Nuc(i)==1.or.Nuc(i)==2.or.Nuc(i)==natom-1.or.Nuc(i)==natom) then
!         end_basis = 1 + end_basis
!      endif
!   end do

!   allocate (iend(2, end_basis)) 
   
   jj=0
   do ii=1, M_in
      if (Nuc(ii)==1.or.Nuc(ii)==natom) then
         jj=jj+1
         Iend(1,jj) = Nuc(ii)
         Iend(2,jj) = ii
      end if
   end do   

end subroutine find_neighbors

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine build_chimera (M_in,fock_in, chimerafock, natom)

   use dftb_data, only:MDFTB, MTB, Iend, end_basis, alfaTB, betaTB, gammaTB

   integer, intent(in)  :: M_in
   integer, intent(in)  :: natom
   real*8, intent(in)   :: fock_in (M_in, M_in)
!#ifdef TD_SIMPLE
   real*8, intent(out)  :: chimerafock (MDFTB, MDFTB) !temporal dimensions
!#else
!   real*16, intent(out) :: chimerafock (2*MTB+M, 2*MTB+M)
!#endif
   integer              :: ii, jj



   chimerafock(:,:) = 0.0D0

   do ii=1, end_basis
      if (Iend(1, ii) == 1) then

         chimerafock(Iend(2,ii)+MTB,MTB) = gammaTB
         chimerafock(MTB, Iend(2,ii)+MTB) =gammaTB
 
      else if (Iend(1, ii) == natom) then

         chimerafock(Iend(2,ii)+MTB,MTB+M_in+1) = gammaTB
         chimerafock(MTB+M_in+1, Iend(2,ii)+MTB) = gammaTB

      end if
   end do


   do ii=1,MTB
      chimerafock(ii,ii) = alfaTB
      chimerafock(MTB+M_in+ii, MTB+M_in+ii)= alfaTB
      
      if (ii<MTB) then

         chimerafock(ii,ii+1) = betaTB
         chimerafock(ii+1,ii) = betaTB 
         chimerafock(MTB+M_in+ii, MTB+M_in+ii+1)= betaTB
         chimerafock(MTB+M_in+ii+1, MTB+M_in+ii)= betaTB

      end if
   end do


   chimerafock(MTB+1:MTB+M_in, MTB+1:MTB+M_in) = fock_in(:,:)


end subroutine build_chimera

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine extract_rhoDFT (M_in, rho_in, rho_out)

   use dftb_data, only: MDFTB, MTB

   implicit none
   integer, intent(in)  :: M_in
   real*8, intent(in)   :: rho_in(MDFTB,MDFTB)
   real*8, intent(out)  :: rho_out(M_in,M_in)
   integer              :: ii, jj

   rho_out=0.0D0

   do jj=1, M_in
   do ii=1, M_in
      rho_out(ii,jj)=rho_in(MTB+ii,MTB+jj)
   end do
   end do

end subroutine extract_rhoDFT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module dftb_subs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
