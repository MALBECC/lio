subroutine truncated_MOs(CoefA,EneA,C_scf,E_scf,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
use excited_data, only: trunc_mos, energy_cut
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
         ! Not truncated MOs.
         call no_trunc(CoefA,EneA,C_scf,E_scf,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
      case ( 1 )
         ! This routine applies the FCA method.
         call fcaApp(CoefA,EneA,C_scf,E_scf,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
      case ( 2 ) 
         ! This routine applies the Reduced Sub-space.
         call reduced_space(CoefA,EneA,C_scf,E_scf,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
      case ( 3 )
         ! This routine delete a specified occupied or virtual MOs.
         call delete_mos(CoefA,EneA,C_scf,E_scf,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
      case default
         print*, "Error in trunc_mos, this only can be 0, 1, 2 or 3"
         stop
   end select

   ! Truncated Range Energy
   ! TODO: for the moment this does not work
   if ( energy_cut ) then
      print*, "WARNING: for the moment this does not work correctly"
      call energy_cutoff(CoefA,EneA,C_scf,E_scf,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
   endif

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

subroutine energy_cutoff(CoefA,EneA,C_scf,E_scf,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
use excited_data,only: map_occ, map_vir, energy_min, energy_max
   implicit none

   integer, intent(in) :: NCO, M
   LIODBLE, intent(in) :: CoefA(:,:), EneA(:)
   integer, intent(inout) :: NCOlr, Mlr, Nvirt, Ndim
   LIODBLE, allocatable, intent(inout) :: C_scf(:,:), E_scf(:)

   integer :: ii, cont, occ, virt, NCOlr_new, Nvirt_new
   LIODBLE :: deltaE
   integer, dimension(:), allocatable :: temp_occ, temp_vir

   print*, "* Range Energy Cutoff" 
   write(*,"(3X,A,F10.5)") "energy_min= ", energy_min
   write(*,"(3X,A,F10.5)") "energy_max= ", energy_max

   ! Get Energy Difference
   allocate(temp_occ(NCOlr),temp_vir(Nvirt))
   temp_occ = 0; temp_vir = 0
   do ii=1,Ndim
      cont = ii - 1
      occ  = NCOlr - (cont/Nvirt)
      virt = mod(cont,Nvirt) + NCOlr + 1
      deltaE = E_scf(virt) - E_scf(occ)
      print*, occ, virt, deltaE
      if ( deltaE > energy_min ) then
      if ( deltaE < energy_max .or. (energy_max < energy_min) ) then
         temp_occ(occ)        = map_occ(occ)
         temp_vir(virt-NCOlr) = map_vir(virt-NCOlr)
      endif
      endif
   enddo
   deallocate(map_occ,map_vir)

   ! Get the new NCOlr and Nvirt
   cont = 0
   do ii=1,NCOlr
      if ( temp_occ(ii) > 0 ) cont = cont + 1
   enddo
   NCOlr_new = cont

   cont = 0
   do ii=1,Nvirt
      if ( temp_vir(ii) > 0 ) cont = cont + 1
   enddo
   Nvirt_new = cont

   ! Get the new maps indexes
   allocate(map_occ(NCOlr_new), map_vir(Nvirt_new))
   cont = 1
   do ii=1,NCOlr
      if ( temp_occ(ii) > 0 ) then
         map_occ(cont) = temp_occ(ii)
         cont = cont + 1
      endif
   enddo
   cont = 1
   do ii=1,Nvirt
      if ( temp_vir(ii) > 0 ) then
         map_vir(cont) = temp_vir(ii)
         cont = cont + 1
      endif
   enddo
   deallocate(temp_occ,temp_vir)

   ! Set the new index
   NCOlr = NCOlr_new; Nvirt = Nvirt_new
   Mlr   = NCOlr + Nvirt
   Ndim  = NCOlr * Nvirt

   ! Get truncated MOs and Energies
   deallocate(C_scf,E_scf)
   allocate(C_scf(M,Mlr),E_scf(Mlr))
   do ii=1,NCOlr
      cont = map_occ(ii)
      C_scf(:,ii) = CoefA(:,cont)
      E_scf(ii)   = EneA(cont)
   enddo
   do ii=1,Nvirt
      cont = NCO + map_vir(ii)
      C_scf(:,NCOlr+ii) = CoefA(:,cont)
      E_scf(NCOlr+ii)   = EneA(cont)
   enddo
end subroutine energy_cutoff

subroutine delete_mos(Cin,Ein,Cout,Eout,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
use excited_data,only: map_occ, map_vir
! In this case the reduced.dat file must be contains 4 lines
! 1) number of occupied MOs to be including
! 2) which occupied MOS, ej: 1 2 3 ...
! 3) number of virtual MOs to be including
! 4) which virtual MOs, ej: 1 2 3 ...
   implicit none

   integer, intent(in) :: NCO, M
   integer, intent(out) :: NCOlr, Mlr, Nvirt, Ndim
   LIODBLE, intent(in) :: Cin(:,:), Ein(:)
   LIODBLE, allocatable, intent(out) :: Cout(:,:), Eout(:)

   logical :: res
   integer :: ii, ind, vir_eff, occ_eff
   integer, dimension(:), allocatable :: mos_occ, mos_vir

   print*, "* Using Deleted MOs"

   ! Reading which MOs are going to be eliminating
   res = .false.
   inquire(file="reduced.dat",exist=res)
   if ( .not. res ) then
      print*, "The file reduced.dat does not exist"
      stop
   endif

   ! Reading Occupied and Virtual MOs
   open(unit=456,file="reduced.dat")
   read(456,*) occ_eff
   allocate(mos_occ(occ_eff))
   read(456,*) (mos_occ(ii), ii=1,occ_eff)
   read(456,*) vir_eff
   allocate(mos_vir(vir_eff))
   read(456,*) (mos_vir(ii), ii=1,vir_eff)
   close(456)

   ! Reports Error
   if ( occ_eff > NCO ) then
      print*, "The number of occ MOs to be including is bigger than &
              & NCO of full system"
      stop
   endif
   if ( vir_eff > (M-NCO) ) then
      print*, "The number of vir MOs to be including is bigger than &
              & NVIRT of full system"
      stop
   endif

   ! Get occupied MOs
   if ( occ_eff == 1 .and. mos_occ(1) == 0 ) then
      print*, "Full NCO to be including"
      NCOlr = NCO
      if(allocated(map_occ)) deallocate(map_occ)
        allocate(map_occ(NCOlr)) ! map_occ(ind small) -> ind big
      do ii=1,NCOlr
         map_occ(ii) = ii
      enddo
   else
      print*, "The occ MOs including are:"
      NCOlr = occ_eff
      if(allocated(map_occ)) deallocate(map_occ)
        allocate(map_occ(NCOlr)) ! map_occ(ind small) -> ind big
      do ii=1,NCOlr
         map_occ(ii) = mos_occ(ii)
         write(*, "(1X,I3)", ADVANCE="NO") mos_occ(ii)
      enddo
      print*, " "
   endif
   deallocate(mos_occ)

   ! Get virtual MOs
   if ( vir_eff == 1 .and. mos_vir(1) == 0 ) then
      print*, "Full NVIRT to be including"
      Nvirt = M - NCO
      if(allocated(map_vir)) deallocate(map_vir)
        allocate(map_vir(Nvirt)) ! map_occ(ind small) -> ind big
      do ii=1,Nvirt
         map_vir(ii) = ii
      enddo
   else
      print*, "The vir MOs including are:"
      Nvirt = vir_eff
      if(allocated(map_vir)) deallocate(map_vir)
        allocate(map_vir(Nvirt)) ! map_vir(ind small) -> ind big
      do ii=1,Nvirt
         map_vir(ii) = mos_vir(ii)
         write(*, "(1X,I3)", ADVANCE="NO") mos_vir(ii)
      enddo
      print*, " "
   endif
   deallocate(mos_vir)

   ! Set the new indexes
   Mlr = NCOlr + Nvirt
   Ndim= NCOlr * Nvirt
   if(allocated(Cout)) deallocate(Cout)
   if(allocated(Eout)) deallocate(Eout)
     allocate(Cout(M,Mlr),Eout(Mlr))

!  Get MOs and Energies truncated
   do ii=1,NCOlr
      ind = map_occ(ii)
      Cout(:,ii) = Cin(:,ind)
      Eout(ii) = Ein(ind)
   enddo
   do ii=1,Nvirt
      ind = NCO + map_vir(ii)
      Cout(:,NCOlr+ii) = Cin(:,ind)
      Eout(NCOlr+ii) = Ein(ind)
   enddo
end subroutine delete_mos
