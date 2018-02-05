subroutine fockbias_oldinit( Natom, nucoforb, sqsmat, fockbias )
   use general_module, only: vector_selection, read_list, atmorb

   implicit none
   integer, intent(in)  :: Natom
   integer, intent(in)  :: nucoforb(:)
   real*8,  intent(in)  :: sqsmat(:,:)
   real*8,  intent(out) :: fockbias(:,:)

   integer, allocatable :: atom_group(:)
   integer, allocatable :: orb_group(:)
   integer, allocatable :: orb_selection(:)

   integer :: Msize
   real*8  :: weight

   Msize = size( sqsmat, 1 )

   if (.not.allocated(atom_group)) then
      allocate(atom_group(Natom))
      call read_list( 'atomgroup', atom_group )
   end if

   if (.not.allocated(orb_group)) then
      allocate(orb_group(Msize))
      call atmorb( atom_group, nucoforb, orb_group )
   end if
   
   if (.not.allocated(orb_selection)) then
      allocate(orb_selection(Msize))
   endif

   fockbias(:,:) = 0.0d0

   weight = 0.195d0
   call vector_selection( 1, orb_group, orb_selection )
   call fterm_biaspot( Msize, sqsmat, orb_selection, weight, fockbias)

   weight=-weight
   call vector_selection( 2, orb_group, orb_selection )
   call fterm_biaspot( Msize, sqsmat, orb_selection, weight, fockbias )

end subroutine fockbias_oldinit
