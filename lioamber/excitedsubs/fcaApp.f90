subroutine fcaApp(Cin,Ein,Cout,Eout,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
! This routine perform Frozen Core and Valence Approximation.
! For the moment, this is broken. TODO

use excited_data, only: FCA, nfo, nfv
   implicit none

   integer, intent(in) :: NCO, M
   integer, intent(out) :: NCOlr, Mlr, Nvirt, Ndim
   double precision, intent(in) :: Cin(:,:), Ein(:)
   double precision, allocatable, intent(out) :: Cout(:,:), Eout(:)

   if (FCA) then
      print*, "FCA doesn't work"
      stop
   else
      print*, "FCA is not used"

      ! SET INDEX
      nfo = 0; nfv = 0
      Nvirt = M - NCO
      NCOlr = NCO
      Mlr = M
      Ndim = Nvirt * NCOlr

      if(.not.allocated(Cout)) allocate(Cout(M,M))
      if(.not.allocated(Eout)) allocate(Eout(M))
      Cout = Cin; Eout = Ein
   endif
end subroutine fcaApp


