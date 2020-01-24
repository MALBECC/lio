subroutine fcaApp(Cin,Ein,Cout,Eout,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
! This routine perform Frozen Core and Valence Approximation.


! INPUTS
! FCA = bool variable to turn on frozen core app.
! nfo = number of frozen occupied orbitals
! nfv = number of frozen virtual orbitals


! For the moment, this routine delete the lowest nfo occupied 
! orbitals and the highest nfv virtuals orbitals

! TODO: selection range orbitals

use excited_data, only: FCA, nfo, nfv
   implicit none

   integer, intent(in) :: NCO, M
   integer, intent(out) :: NCOlr, Mlr, Nvirt, Ndim
   double precision, intent(in) :: Cin(:,:), Ein(:)
   double precision, allocatable, intent(out) :: Cout(:,:), Eout(:)

   integer :: ii

   if (FCA) then

      print*, "*Using Frozen Core Approximation"
      print*, "number of frozen occupied orbitals", nfo
      print*, "number of frozen virtual orbitals", nfv
      NCOlr = NCO - nfo
      Nvirt = M - NCO - nfv
      Mlr   = M - nfo - nfv
      Ndim  = Nvirt * NCOlr

      if(.not.allocated(Cout)) allocate(Cout(M,Mlr))
      if(.not.allocated(Eout)) allocate(Eout(Mlr))

      ! Set occupied
      do ii=1,NCOlr
         Cout(:,NCOlr+1-ii) = Cin(:,NCO+1-ii)
         Eout(NCOlr+1-ii)   = Ein(NCO+1-ii)
      enddo

      ! Set virtual
      do ii=1,Nvirt
         Cout(:,NCOlr+ii) = Cin(:,NCO+ii)
         Eout(NCOlr+ii)   = Ein(NCO+ii)
      enddo

   else
      print*, "*FCA is not used"

      ! SET INDEX
      nfo   = 0; nfv = 0
      Nvirt = M - NCO
      NCOlr = NCO
      Mlr   = M
      Ndim  = Nvirt * NCOlr

      if(.not.allocated(Cout)) allocate(Cout(M,M))
      if(.not.allocated(Eout)) allocate(Eout(M))

      Cout = Cin; Eout = Ein
   endif
end subroutine fcaApp


