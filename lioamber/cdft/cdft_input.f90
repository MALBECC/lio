! Reads a CDFT input from input file.
subroutine cdft_input_read(input_UID)
   use cdft_data, only: cdft_spin, cdft_chrg, cdft_reg
   implicit none
   integer, intent(in) :: input_UID
   
   character(len=10) :: buffer
   integer           :: ios, inp_spin, inp_chrg, ii
   
   rewind(input_UID)
   ios = 0
   do while ((trim(buffer) /= "{CDFT}") .and. (ios == 0) )
      read(input_UID,'(A10)', iostat=ios) buffer
   enddo

   ! If ios < 0, found EOF. No CDFT input provided.
   if (ios < 0) return

   write(*,'(A)') "CDFT input found, reading options."
   ! Starts reading CDFT data.
   read(input_UID,*) cdft_reg%n_regions, inp_chrg, inp_spin
   if (inp_chrg == 1) cdft_chrg = .true.
   if (inp_spin == 1) cdft_spin = .true.
   write(*,'(A21,I3,A21,L2,A19,L2)')"  Number of regions: ",cdft_reg%n_regions,&
                                    " | Constrain charge: ", cdft_chrg, &
                                    " | Constrain spin: ", cdft_spin

   if (allocated(cdft_reg%chrg))  deallocate(cdft_reg%chrg)
   if (allocated(cdft_reg%spin))  deallocate(cdft_reg%spin)
   if (allocated(cdft_reg%natom)) deallocate(cdft_reg%natom)
   allocate(cdft_reg%chrg(cdft_reg%n_regions))
   allocate(cdft_reg%spin(cdft_reg%n_regions))
   allocate(cdft_reg%natom(cdft_reg%n_regions))

   do ii = 1, cdft_reg%n_regions
      read(input_UID,*) cdft_reg%natom(ii), cdft_reg%chrg(ii), &
                        cdft_reg%spin(ii)
   enddo
   if (allocated(cdft_reg%atoms)) deallocate(cdft_reg%atoms)
   cdft_reg%max_nat = maxval(cdft_reg%natom,1)
   allocate(cdft_reg%atoms(cdft_reg%n_regions, cdft_reg%max_nat))
   cdft_reg%atoms = 0
   do ii = 1, cdft_reg%n_regions
      read(input_UID,*) cdft_reg%atoms(ii,1:cdft_reg%natom(ii))
   enddo
end subroutine cdft_input_read

! Checks if doing CDFT and sets energy_all_iterations to true in order to
! also calculate Becke charges in each iteration step.
subroutine cdft_options_check(do_becke, open_shell)
   use cdft_data, only: cdft_chrg, cdft_spin, doing_cdft

   implicit none
   logical, intent(inout) :: do_becke, open_shell

   if ((cdft_chrg) .or. (cdft_spin)) then
      doing_cdft     = .true.
      do_becke       = .true.
   endif

   if (cdft_spin) open_shell = .true.
end subroutine cdft_options_check