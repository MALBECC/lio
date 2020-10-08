! Initializes arrays.
subroutine cdft_initialise(n_atoms)
   use cdft_data, only: cdft_reg, cdft_c

   implicit none
   integer, intent(in) :: n_atoms
   integer             :: J_size = 0

   if (cdft_c%do_chrg .and. cdft_c%do_spin) then
      J_size = 2 * cdft_c%n_regions
      cdft_c%sp_idx = cdft_c%n_regions
   else if (cdft_c%do_chrg .or. cdft_c%do_spin) then
      J_size = cdft_c%n_regions
      cdft_c%sp_idx = 0
   endif

   if (cdft_c%do_chrg .or. cdft_c%do_spin) then
      if (allocated(cdft_c%at_chrg))  deallocate(cdft_c%at_chrg)
      if (allocated(cdft_reg%Vc))     deallocate(cdft_reg%Vc)
      if (allocated(cdft_reg%Vc_old)) deallocate(cdft_reg%Vc_old)
      allocate(cdft_c%at_chrg(n_atoms))
      allocate(cdft_reg%Vc(cdft_c%n_regions))
      allocate(cdft_reg%Vc_old(cdft_c%n_regions))

      if (allocated(cdft_c%at_spin))  deallocate(cdft_c%at_spin)
      if (allocated(cdft_reg%Vs))     deallocate(cdft_reg%Vs)
      if (allocated(cdft_reg%Vs_old)) deallocate(cdft_reg%Vs_old)
      allocate(cdft_c%at_spin(n_atoms))  
      allocate(cdft_reg%Vs(cdft_c%n_regions))
      allocate(cdft_reg%Vs_old(cdft_c%n_regions))

      if (allocated(cdft_reg%cst))     deallocate(cdft_reg%cst)
      if (allocated(cdft_reg%cst_old)) deallocate(cdft_reg%cst_old)
      if (allocated(cdft_reg%Vmix))    deallocate(cdft_reg%Vmix)
      if (allocated(cdft_reg%Vm_old))  deallocate(cdft_reg%Vm_old)
      allocate(cdft_reg%cst(J_size))
      allocate(cdft_reg%cst_old(J_size))
      allocate(cdft_reg%Vmix(J_size))
      allocate(cdft_reg%Vm_old(J_size))

      if (cdft_c%mixed) then
         if (allocated(cdft_reg%Vc2))  deallocate(cdft_reg%Vc2)
         if (allocated(cdft_reg%Vs2))  deallocate(cdft_reg%Vs2)
         allocate(cdft_reg%Vc2(cdft_c%n_regions))
         allocate(cdft_reg%Vs2(cdft_c%n_regions))
         cdft_reg%Vc2 = 0.0D0                      
         cdft_reg%Vs2 = 0.0D0      
      endif
   endif

   if (allocated(cdft_c%jacob)) deallocate(cdft_c%jacob)
   allocate(cdft_c%jacob(J_size, J_size))

   call g2g_cdft_init(cdft_c%do_chrg, cdft_c%do_spin, cdft_c%n_regions, &
                      cdft_c%max_nat, cdft_reg%natom, cdft_reg%atoms)
   cdft_reg%Vc = 0.0D0                      
   cdft_reg%Vs = 0.0D0
   
   call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
end subroutine cdft_initialise

! Deallocates arrays.
subroutine cdft_finalise()
   use cdft_data, only: cdft_reg, cdft_c
   implicit none

   cdft_c%do_chrg   = .false.
   cdft_c%do_spin   = .false.
   cdft_c%sp_idx    = 0
   cdft_c%max_nat   = 0
   cdft_c%n_regions = 0

   if (allocated(cdft_c%jacob))     deallocate(cdft_c%jacob)
   if (allocated(cdft_c%at_chrg))   deallocate(cdft_c%at_chrg)
   if (allocated(cdft_c%at_spin))   deallocate(cdft_c%at_spin)
   if (allocated(cdft_reg%Vc))      deallocate(cdft_reg%Vc)
   if (allocated(cdft_reg%Vs))      deallocate(cdft_reg%Vs)
   if (allocated(cdft_reg%Vc2))     deallocate(cdft_reg%Vc2)
   if (allocated(cdft_reg%Vs2))     deallocate(cdft_reg%Vs2)
   if (allocated(cdft_reg%natom))   deallocate(cdft_reg%natom)
   if (allocated(cdft_reg%atoms))   deallocate(cdft_reg%atoms)
   if (allocated(cdft_reg%chrg))    deallocate(cdft_reg%chrg)
   if (allocated(cdft_reg%spin))    deallocate(cdft_reg%spin)
   if (allocated(cdft_reg%chrg2))   deallocate(cdft_reg%chrg2)
   if (allocated(cdft_reg%spin2))   deallocate(cdft_reg%spin2)
   if (allocated(cdft_reg%Vc_old))  deallocate(cdft_reg%Vc_old)
   if (allocated(cdft_reg%Vs_old))  deallocate(cdft_reg%Vs_old)
   if (allocated(cdft_reg%Vmix))    deallocate(cdft_reg%Vmix)
   if (allocated(cdft_reg%Vm_old))  deallocate(cdft_reg%Vm_old)
   if (allocated(cdft_reg%cst))     deallocate(cdft_reg%cst)
   if (allocated(cdft_reg%cst_old)) deallocate(cdft_reg%cst_old)

   call g2g_cdft_finalise()
end subroutine cdft_finalise
