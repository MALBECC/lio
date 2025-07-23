subroutine fcaApp(Cin,Ein,Cout,Eout,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
! This routine perform Frozen Core and Valence Approximation.


! INPUTS
! FCA = bool variable to turn on frozen core app.
! nfo = number of frozen occupied orbitals
! nfv = number of frozen virtual orbitals


! For the moment, this routine delete the lowest nfo occupied 
! orbitals and the highest nfv virtuals orbitals

! TODO: selection range orbitals

use excited_data, only: nfo, nfv, map_occ, map_vir

   implicit none

   integer, intent(in) :: NCO, M
   integer, intent(out) :: NCOlr, Mlr, Nvirt, Ndim
   LIODBLE, intent(in) :: Cin(:,:), Ein(:)
   LIODBLE, allocatable, intent(out) :: Cout(:,:), Eout(:)

   integer :: ii

   print*, "*Using Frozen Core Approximation"
   print*, "number of frozen occupied orbitals", nfo
   print*, "number of frozen virtual orbitals", nfv

   ! Set indexes
   NCOlr = NCO - nfo
   Nvirt = M - NCO - nfv
   Mlr   = M - nfo - nfv
   Ndim  = Nvirt * NCOlr

   if(.not.allocated(Cout)) allocate(Cout(M,Mlr))
   if(.not.allocated(Eout)) allocate(Eout(Mlr))
   if(allocated(map_occ)) deallocate(map_occ)
   if(allocated(map_vir)) deallocate(map_vir)
     allocate(map_occ(NCOlr),map_vir(Nvirt))

   ! Set occupied
   do ii=1,NCOlr
      Cout(:,NCOlr+1-ii) = Cin(:,NCO+1-ii)
      Eout(NCOlr+1-ii)   = Ein(NCO+1-ii)
      map_occ(NCOlr+1-ii)= NCO+1-ii
   enddo

   ! Set virtual
   do ii=1,Nvirt
      Cout(:,NCOlr+ii) = Cin(:,NCO+ii)
      Eout(NCOlr+ii)   = Ein(NCO+ii)
      map_vir(ii)      = ii
   enddo
end subroutine fcaApp

subroutine fca_restored(Cin,Ein,Cout,Eout,Xin,M,Mlr,Nvirt,NCO,NCOlr,Ndim,Nstat)
use excited_data, only: nfo, nfv, map_occ, map_vir, trunc_mos
   implicit none

   integer, intent(in)  :: M, NCO, Nstat
   LIODBLE, intent(in) :: Cin(:,:), Ein(:)

   integer, intent(inout) :: NCOlr, Mlr, Nvirt, Ndim
   LIODBLE, allocatable, intent(inout) :: Cout(:,:), Eout(:), Xin(:,:)

   integer :: ii, jj, kk, ind1, ind2, ii_fix, jj_fix
   integer :: Nvirt_out, Ndim_out
   LIODBLE, allocatable :: Xout(:,:)

   if ( trunc_mos == 0 ) return 
   print*, "Restoring Orbitals"

   ! Deinitialization
   call basis_deinitLR()

   ! Set new parameter
   Nvirt_out = M - NCO
   Ndim_out  = NCO * Nvirt_out
  
   ! Form change basis matrix with full MO
   call basis_initLR(Cin,M,M,NCO,Nvirt_out)
   deallocate(Cout,Eout); allocate(Cout(M,M),Eout(M))
   Cout = Cin; Eout = Ein

   ! New Transition Vectors with FULL DIMENSION
   allocate(Xout(Ndim_out,Nstat)); Xout = 0.0d0
   do ii=1,NCOlr
   do jj=1,Nvirt
      ii_fix = NCO - map_occ(NCOlr+1-ii) + 1
      jj_fix = map_vir(jj)
      ind1 = (ii-1) * Nvirt + jj
      ind2 = (ii_fix-1) * Nvirt_out + jj_fix

      do kk=1,Nstat
         Xout(ind2,kk) = Xin(ind1,kk)
      enddo
   enddo
   enddo

   deallocate(Xin); allocate(Xin(Ndim_out,Nstat))
   Xin = Xout
   deallocate(Xout)

   ! Final indexes
   NCOlr = NCO; Mlr = M; Nvirt = Nvirt_out; Ndim = Ndim_out
   nfo = 0; nfv = 0
end subroutine fca_restored
