! Reads a CDFT input from input file.
subroutine cdft_input_read(input_UID)
   use cdft_data, only: cdft_c, cdft_reg
   implicit none
   integer, intent(in) :: input_UID
   
   character(len=10) :: buffer
   integer           :: ios = 0
   integer           :: inp_spin  = 0
   integer           :: inp_chrg  = 0
   integer           :: inp_mixed = 0
   integer           :: ii
   
   rewind(input_UID)
   ios = 0
   do while ((trim(buffer) /= "{CDFT}") .and. (ios == 0) )
      read(input_UID,'(A10)', iostat=ios) buffer
   enddo

   ! If ios < 0, found EOF. No CDFT input provided.
   if (ios < 0) return

   write(*,'(A)') "CDFT input found, reading options."
   ! Starts reading CDFT data. It tries to read four options,
   ! accounting for CDFT; if not, it reads three for a regular
   ! CDFT calculation.
   ios = 0
   read(input_UID,*,iostat=ios) cdft_c%n_regions, inp_chrg, &
                                inp_spin, inp_mixed
   if (ios /= 0) read(input_UID,*) cdft_c%n_regions, inp_chrg, inp_spin

   if (inp_mixed == 1) cdft_c%mixed   = .true.
   if (inp_chrg  == 1) cdft_c%do_chrg = .true.
   if (inp_spin  == 1) cdft_c%do_spin = .true.
   write(*,'(A21,I3,A21,L2,A19,L2)')"  Number of regions: ",cdft_c%n_regions,&
                                    " | Constrain charge: ", cdft_c%do_chrg, &
                                    " | Constrain spin: "  , cdft_c%do_spin

   if (allocated(cdft_reg%chrg))  deallocate(cdft_reg%chrg)
   if (allocated(cdft_reg%spin))  deallocate(cdft_reg%spin)
   if (allocated(cdft_reg%natom)) deallocate(cdft_reg%natom)
   allocate(cdft_reg%chrg(cdft_c%n_regions))
   allocate(cdft_reg%spin(cdft_c%n_regions))
   allocate(cdft_reg%natom(cdft_c%n_regions))

   if (cdft_c%mixed) then
      write(*,'(A)') "  Doing mixed CDFT for Hab calculation."
      if (allocated(cdft_reg%chrg2)) deallocate(cdft_reg%chrg2)
      if (allocated(cdft_reg%spin2)) deallocate(cdft_reg%spin2)
      allocate(cdft_reg%chrg2(cdft_c%n_regions))
      allocate(cdft_reg%spin2(cdft_c%n_regions))
   endif

   if (cdft_c%mixed) then
      do ii = 1, cdft_c%n_regions
         read(input_UID,*) cdft_reg%natom(ii), cdft_reg%chrg(ii), &
                           cdft_reg%chrg2(ii), cdft_reg%spin(ii), &
                           cdft_reg%spin2(ii)
      enddo
   else
      do ii = 1, cdft_c%n_regions
        read(input_UID,*) cdft_reg%natom(ii), cdft_reg%chrg(ii), &
                          cdft_reg%spin(ii)
      enddo
   endif
   if (allocated(cdft_reg%atoms)) deallocate(cdft_reg%atoms)
   cdft_c%max_nat = maxval(cdft_reg%natom,1)
   allocate(cdft_reg%atoms(cdft_c%n_regions, cdft_c%max_nat))
   cdft_reg%atoms = 0
   do ii = 1, cdft_c%n_regions
      read(input_UID,*) cdft_reg%atoms(ii,1:cdft_reg%natom(ii))
   enddo

   rewind(input_UID)
end subroutine cdft_input_read

! Checks if doing CDFT and sets energy_all_iterations to true in order to
! also calculate Becke charges in each iteration step.
subroutine cdft_options_check(do_becke, open_shell)
   use cdft_data, only: cdft_c, doing_cdft

   implicit none
   logical, intent(inout) :: do_becke, open_shell

   if ((cdft_c%do_chrg) .or. (cdft_c%do_spin)) then
      doing_cdft = .true.
      do_becke   = .true.
   endif

   if (cdft_c%do_spin) open_shell = .true.
end subroutine cdft_options_check