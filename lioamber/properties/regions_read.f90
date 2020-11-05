! Reads region inputs, in order to perform region printing.
! The region input should have the following structure:
! {PROPREGIONS}
!  NREGS
!  REG1_NATOMS REG2_NATOMS REGN_NATOMS
!  REG1_ATOM1  REG1_ATOM2  REG1_ATOMN
!  REG2_ATOM1  REG2_ATOM2  REG2_ATOMN
!  REGN_ATOM1  REGN_ATOM2  REGN_ATOMN
! {END}
subroutine properties_region_read(input_UID)
   use properties_data, only: prop_regions
   implicit none
   integer, intent(in) :: input_UID
   
   character(len=20) :: buffer
   integer           :: ios, ii, max_nat
   
   rewind(input_UID)
   ios = 0
   buffer = " "
   do while ((trim(buffer) /= "{PROPREGIONS}") .and. (ios == 0) )
      read(input_UID,'(A20)', iostat=ios) buffer
   enddo

   ! If ios < 0, found EOF. No region input provided.
   if (ios < 0) return

   ! Starts reading region data.
   read(input_UID,*) prop_regions%n_regions
   allocate(prop_regions%natoms(prop_regions%n_regions))

   read(input_UID,*) prop_regions%natoms(1:prop_regions%n_regions)

   max_nat = maxval(prop_regions%natoms,1)
   if (allocated(prop_regions%atoms)) deallocate(prop_regions%atoms)
   allocate(prop_regions%atoms(prop_regions%n_regions, max_nat))
   prop_regions%atoms = 0
   do ii = 1, prop_regions%n_regions
      read(input_UID,*) prop_regions%atoms(ii,1:prop_regions%natoms(ii))
   enddo
   rewind(input_UID)
end subroutine properties_region_read
