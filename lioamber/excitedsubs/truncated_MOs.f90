subroutine truncated_MOs(CoefA,EneA,C_scf,E_scf,map_occ,map_vir,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
use excited_data, only: trunc_mos
use garcha_mod  , only: OPEN

! This routine perform truncated Molecular Orbitals with FCA or Reduced Subspace
! VARIABLES
   ! NCO = number of occupied molecular orbitals
   ! Nvirt = number of virtual molecular orbitals
   ! Ndim = dimension of Excited Matrix = (NCOxNvirt)^2
   ! C_scf, E_scf = Molecular Orbital Coeff. and Energy
   ! map_occ = map the actual index in C_scf to the real MO occupied in CoefA
   ! map_vir = map the actual index in C_scf to the real MO virtual in CoefA
   implicit none

   integer, intent(in)  :: NCO, M
   integer, intent(out) :: NCOlr, Mlr, Nvirt, Ndim
   LIODBLE, intent(in)  :: CoefA(:,:), EneA(:)
   integer, allocatable, intent(out) :: map_occ(:), map_vir(:)
   LIODBLE, allocatable, intent(out) :: C_scf(:,:), E_scf(:)

   if ( OPEN .and. trunc_mos /= 0 ) then
      print*, "Linear Response Open shell only works with trunc_mos=0 option"
      stop
   endif

   ! Truncated Molecular Orbitals
   select case ( trunc_mos ) 
      case ( 0 )
         ! Not truncated MOs.
         call no_trunc(CoefA,EneA,C_scf,E_scf,map_occ,map_vir,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
      case ( 1 )
         ! This routine applies the FCA method.
         call fcaApp(CoefA,EneA,C_scf,E_scf,map_occ,map_vir,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
      case ( 2 ) 
         ! This routine applies the Atoms Reduced Sub-space.
         call reduced_space(CoefA,EneA,C_scf,E_scf,map_occ,map_vir,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
      case ( 3 )
         ! This routine delete a specified occupied or virtual MOs.
         call delete_mos(CoefA,EneA,C_scf,E_scf,map_occ,map_vir,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
      case default
         print*, "Error in trunc_mos, this only can be 0, 1, 2 or 3"
         stop
   end select

end subroutine truncated_MOs

subroutine no_trunc(Cin,Ein,Cout,Eout,map_occ,map_vir,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
   implicit none

   integer, intent(in) :: NCO, M
   integer, intent(out) :: NCOlr, Mlr, Nvirt, Ndim
   LIODBLE, intent(in) :: Cin(:,:), Ein(:)
   integer, allocatable, intent(out) :: map_occ(:), map_vir(:)
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

subroutine delete_mos(Cin,Ein,Cout,Eout,map_occ,map_vir,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
! In this case the reduced.dat file must contain 4 lines
! 1) number of occupied MOs to be including
! 2) which occupied MOS, ej: 1 2 3 ...
! 3) number of virtual MOs to be including
! 4) which virtual MOs, ej: 1 2 3 ...
   implicit none

   integer, intent(in)  :: NCO, M
   integer, intent(out) :: NCOlr, Mlr, Nvirt, Ndim
   LIODBLE, intent(in)  :: Cin(:,:), Ein(:)
   integer, allocatable, intent(out) :: map_occ(:), map_vir(:)
   LIODBLE, allocatable, intent(out) :: Cout(:,:), Eout(:)

   logical :: res
   integer :: ii, ind, vir_eff, occ_eff
   integer, dimension(:), allocatable :: mos_occ, mos_vir

   print*, "*Using Deleted MOs"

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
         write(*, "(1X,I3)", ADVANCE="NO") NCO+mos_vir(ii)
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
