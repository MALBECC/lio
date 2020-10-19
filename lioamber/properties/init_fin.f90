subroutine properties_initialise(open_shell, do_td, atom_of_func)
   use properties_data, only: dipole, prop_regions
   implicit none

   logical, intent(in) :: open_shell
   integer, intent(in) :: do_td
   integer, intent(in) :: atom_of_func(:)

   call properties_open_all(open_shell, do_td, 0)

   if (prop_regions%n_regions > 0) then
      call properties_open_all(open_shell, do_td, 50)
      if (dipole .or. (do_td > 0)) &
            call initialise_regions_functions(atom_of_func)
   endif
end subroutine properties_initialise

! Opens all files, using an offset for separate-region output.
subroutine properties_open_all(open_shell, do_td, offset)
   use fileio         , only: safeio_open
   use properties_data, only: fmulliken, fdipole, flowdin, fbecke,     &
                              ffukui, mulliken, dipole, lowdin, becke, &
                              fukui, UIDs
   implicit none
   logical, intent(in) :: open_shell
   integer, intent(in) :: do_td
   integer, intent(in) :: offset
   character(len=60)   :: tmp_name
   character(len=7)    :: region
   
   region = "region_"

   if (offset == 0) then
      if (mulliken) then
         call safeio_open(UIDs%mul, fmulliken, 3)
         tmp_name = trim(fmulliken)//"_spin"
         if (open_shell) call safeio_open(UIDs%muls, tmp_name, 3)
      endif

      if (lowdin) then
         call safeio_open(UIDs%low, flowdin, 3)
         tmp_name = trim(flowdin)//"_spin"
         if (open_shell) call safeio_open(UIDs%lows, tmp_name, 3)
      endif

      if (becke) then
         call safeio_open(UIDs%bec, fbecke, 3)
         tmp_name = trim(fbecke)//"_spin"
         if (open_shell) call safeio_open(UIDs%becs, tmp_name, 3)
      endif

      if (fukui) call safeio_open(UIDs%fuk, ffukui, 3)
      if (dipole) call safeio_open(UIDs%dip, fdipole, 3)

      if (do_td > 0) then
         tmp_name = trim(fdipole)//"_td"
         call safeio_open(UIDs%diptd, tmp_name, 3)
      endif
   else
      if (mulliken) then
         call safeio_open(UIDs%mul+offset, trim(region)//fmulliken, 3)
         tmp_name = trim(region)//trim(fmulliken)//"_spin"
         if (open_shell) call safeio_open(UIDs%muls+offset, tmp_name, 3)
      endif

      if (lowdin) then
         call safeio_open((UIDs%low+offset), trim(region)//flowdin, 3)
         tmp_name = trim(region)//trim(flowdin)//"_spin"
         if (open_shell) call safeio_open(UIDs%lows+offset, tmp_name, 3)
      endif

      if (becke) then
         call safeio_open((UIDs%bec+offset), trim(region)//fbecke, 3)
         tmp_name = trim(region)//trim(fbecke)//"_spin"
         if (open_shell) call safeio_open(UIDs%becs+offset, tmp_name, 3)
      endif

      if (dipole) call safeio_open((UIDs%dip+offset), trim(region)//fdipole, 3)

      if (do_td > 0) then
         tmp_name = trim(region)//trim(fdipole)//"_td"
         call safeio_open(UIDs%diptd+offset, tmp_name, 3)
      endif
   endif

end subroutine properties_open_all

! Reorders data to get which functions correspond to each region.
! Funcion indexes are stored in %funcs(region_id, n_function)
subroutine initialise_regions_functions(atom_of_func)
   use properties_data, only: prop_regions
   implicit none

   integer, intent(in) :: atom_of_func(:)
   integer :: ifunct, ireg, max_funcs, f_index
   integer :: vec_size


   allocate(prop_regions%nfuncs(prop_regions%n_regions))
   prop_regions%nfuncs = 0

   do ifunct = 1, size(atom_of_func,1)
      do iReg = 1, prop_regions%n_regions
         if (any(prop_regions%atoms(iReg,:) == atom_of_func(ifunct))) &
               prop_regions%nfuncs(iReg) = prop_regions%nfuncs(iReg) +1 
      enddo
   enddo

   max_funcs = maxval(prop_regions%nfuncs,1)
   allocate(prop_regions%funcs(prop_regions%n_regions, max_funcs))
   prop_regions%funcs = 0

   do iReg = 1, prop_regions%n_regions
      f_index = 1

      do ifunct = 1, size(atom_of_func)
         if (any(prop_regions%atoms(iReg,:) == atom_of_func(ifunct))) then
            prop_regions%funcs(iReg,f_index) = ifunct
            f_index = f_index +1 
         endif
      enddo
   enddo
end subroutine initialise_regions_functions

subroutine properties_finalise()
   use properties_data, only: prop_regions
   implicit none

   ! Closes all regular output files.
   call properties_close_all(0)

   ! Closes all region-related output files.
   call properties_close_all(50)

   prop_regions%n_regions = 0
   if (allocated(prop_regions%natoms)) deallocate(prop_regions%natoms)
   if (allocated(prop_regions%atoms) ) deallocate(prop_regions%atoms)
   if (allocated(prop_regions%nfuncs)) deallocate(prop_regions%nfuncs)
   if (allocated(prop_regions%funcs) ) deallocate(prop_regions%funcs)

end subroutine properties_finalise

subroutine properties_close_all(offset)
   use fileio         , only: safeio_close
   use properties_data, only: UIDs
   implicit none
   integer, intent(in) :: offset

   call safeio_close(UIDs%dip   + offset)
   call safeio_close(UIDs%diptd + offset)
   call safeio_close(UIDs%fuk   + offset)
   call safeio_close(UIDs%bec   + offset)
   call safeio_close(UIDs%becs  + offset)
   call safeio_close(UIDs%low   + offset)
   call safeio_close(UIDs%lows  + offset)
   call safeio_close(UIDs%mul   + offset)
   call safeio_close(UIDs%muls  + offset)

end subroutine properties_close_all