! Initializes arrays.
subroutine cdft_initialise(n_atoms)
   use cdft_data, only: cdft_reg, cdft_chrg, cdft_spin, &
                        at_chrg, at_spin, jacob, sp_idx
   implicit none
   integer, intent(in) :: n_atoms
   integer             :: J_size = 0

   if (cdft_chrg .and. cdft_spin) then
      J_size = 2 * cdft_reg%n_regions
      sp_idx = cdft_reg%n_regions
   else if (cdft_chrg .or. cdft_spin) then
      J_size = cdft_reg%n_regions
      sp_idx = 0
   endif

   if (cdft_chrg .or. cdft_spin) then
      if (allocated(at_chrg))         deallocate(at_chrg)
      if (allocated(cdft_reg%Vc))     deallocate(cdft_reg%Vc)
      if (allocated(cdft_reg%Vc_old)) deallocate(cdft_reg%Vc_old)
      allocate(at_chrg(n_atoms))
      allocate(cdft_reg%Vc(cdft_reg%n_regions))
      allocate(cdft_reg%Vc_old(cdft_reg%n_regions))

      if (allocated(at_spin))         deallocate(at_spin)
      if (allocated(cdft_reg%Vs))     deallocate(cdft_reg%Vs)
      if (allocated(cdft_reg%Vs_old)) deallocate(cdft_reg%Vs_old)
      allocate(at_spin(n_atoms))  
      allocate(cdft_reg%Vs(cdft_reg%n_regions))
      allocate(cdft_reg%Vs_old(cdft_reg%n_regions))

      if (allocated(cdft_reg%cst))     deallocate(cdft_reg%cst)
      if (allocated(cdft_reg%cst_old)) deallocate(cdft_reg%cst_old)
      if (allocated(cdft_reg%Vmix))    deallocate(cdft_reg%Vmix)
      if (allocated(cdft_reg%Vm_old))  deallocate(cdft_reg%Vm_old)
      allocate(cdft_reg%cst(J_size))
      allocate(cdft_reg%cst_old(J_size))
      allocate(cdft_reg%Vmix(J_size))
      allocate(cdft_reg%Vm_old(J_size))
   endif

   if (allocated(jacob)) deallocate(jacob)
   allocate(jacob(J_size, J_size))

   call g2g_cdft_init(cdft_chrg, cdft_spin, cdft_reg%n_regions, &
                      cdft_reg%max_nat, cdft_reg%natom, cdft_reg%atoms)
   cdft_reg%Vc = 0.0D0                      
   cdft_reg%Vs = 0.0D0
   call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
end subroutine cdft_initialise

! Deallocates arrays.
subroutine cdft_finalise()
   use cdft_data, only: cdft_reg, at_spin, at_chrg, jacob
   implicit none

   if (allocated(jacob))            deallocate(jacob)
   if (allocated(at_chrg))          deallocate(at_chrg)
   if (allocated(at_spin))          deallocate(at_spin)
   if (allocated(cdft_reg%Vc))      deallocate(cdft_reg%Vc)
   if (allocated(cdft_reg%Vs))      deallocate(cdft_reg%Vs)
   if (allocated(cdft_reg%natom))   deallocate(cdft_reg%natom)
   if (allocated(cdft_reg%atoms))   deallocate(cdft_reg%atoms)
   if (allocated(cdft_reg%chrg))    deallocate(cdft_reg%chrg)
   if (allocated(cdft_reg%spin))    deallocate(cdft_reg%spin)
   if (allocated(cdft_reg%Vc_old))  deallocate(cdft_reg%Vc_old)
   if (allocated(cdft_reg%Vs_old))  deallocate(cdft_reg%Vs_old)
   if (allocated(cdft_reg%Vmix))    deallocate(cdft_reg%Vmix)
   if (allocated(cdft_reg%Vm_old))  deallocate(cdft_reg%Vm_old)
   if (allocated(cdft_reg%cst))     deallocate(cdft_reg%cst)
   if (allocated(cdft_reg%cst_old)) deallocate(cdft_reg%cst_old)

   call g2g_cdft_finalise()
end subroutine cdft_finalise