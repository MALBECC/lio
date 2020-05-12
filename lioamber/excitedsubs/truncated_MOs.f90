subroutine truncated_MOs(CoefA,EneA,C_scf,E_scf,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
use excited_data, only: trunc_mos
! This routine perform truncated Molecular Orbitals with FCA or Reduced Subspace
! VARIABLES
   ! NCO = number of occupied molecular orbitals
   ! Nvirt = number of virtual molecular orbitals
   ! Ndim = dimension of Excited Matrix = (NCOxNvirt)^2
   ! C_scf, E_scf = Molecular Orbital Coeff. and Energy
   implicit none

   integer, intent(in)  :: NCO, M
   integer, intent(out) :: NCOlr, Mlr, Nvirt, Ndim
   LIODBLE, intent(in)  :: CoefA(:,:), EneA(:)
   LIODBLE, allocatable, intent(out) :: C_scf(:,:), E_scf(:)

   ! Truncated Molecular Orbitals
   select case ( trunc_mos ) 
      case ( 0 )
         ! Not truncated MOs
         call no_trunc(CoefA,EneA,C_scf,E_scf,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
      case ( 1 )
         ! This routine applies the FCA method
         call fcaApp(CoefA,EneA,C_scf,E_scf,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
      case ( 2 ) 
         ! This routine applies the Reduced Sub-space.
         call reduced_space(CoefA,EneA,C_scf,E_scf,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
      case default
         print*, "Error in trunc_mos, this only can be 0, 1 or 2"
         stop
   end select

   ! TODO: Truncated Range Energy

end subroutine truncated_MOs

subroutine no_trunc(Cin,Ein,Cout,Eout,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
use excited_data, only: map_occ, map_vir
   implicit none

   integer, intent(in) :: NCO, M
   integer, intent(out) :: NCOlr, Mlr, Nvirt, Ndim
   LIODBLE, intent(in) :: Cin(:,:), Ein(:)
   LIODBLE, allocatable, intent(out) :: Cout(:,:), Eout(:)

   integer :: ii
   print*, "*Truncated MOs is not used"

   ! SET INDEX
   Nvirt = M - NCO
   NCOlr = NCO
   Mlr   = M
   Ndim  = Nvirt * NCOlr

   if(.not.allocated(Cout)) allocate(Cout(M,M))
   if(.not.allocated(Eout)) allocate(Eout(M))
   Cout = Cin; Eout = Ein

   ! Get small index to big index
   if(allocated(map_occ)) deallocate(map_occ)
   if(allocated(map_vir)) deallocate(map_vir)
     allocate(map_occ(NCO),map_vir(Nvirt))
   do ii=1,NCO
      map_occ(ii) = ii
   enddo
   do ii=1,Nvirt
      map_vir(ii) = ii
   enddo

end subroutine no_trunc
