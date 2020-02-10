!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Check_posvel( time_step, natoms, posvec, velvec, fname )
   implicit none
   real*8          , intent(in) :: time_step
   integer         , intent(in) :: natoms
   real*8          , intent(in) :: posvec(3, natoms)
   real*8          , intent(in) :: velvec(3, natoms)
   character(len=*), intent(in) :: fname

   integer                     :: funit
   integer                     :: int_stat
   real*8                      :: posdif
   real*8                      :: veltry
   integer                     :: natom, kdir
   real*8, allocatable, save   :: posold(:,:)
   logical            , save   :: first_call = .true.
   character(len=*), parameter :: strfmt = "(2(2x,I6),5(2x,E15.8))"

#  ifndef DEBUGGING
      return
#  endif
   if (.not.allocated(posold)) then
      allocate( posold(3,natoms) )
      posold(:,:) = 0.0d0
   endif

   funit = 1234
   if (first_call) then
      first_call = .false.
      open( unit=funit, file=fname, iostat=int_stat )
   else
      open( unit=funit, file=fname, position="append", iostat=int_stat )
   endif

   if ( int_stat /= 0 ) then
      print*, "There is something wrong with output file..."
      print*, "file = ", fname
      print*, "unit = ", funit
      print*, "iost = ", int_stat
      print*, ""
   end if


   do natom = 1, natoms
   do kdir = 1, 3
      posdif = posvec(kdir, natom) - posold(kdir, natom)
      veltry = posdif / time_step
      if ( int_stat == 0 ) then
         write( unit=funit, fmt=strfmt ) natom, kdir, &
         & posold(kdir, natom), posvec(kdir, natom), posdif, &
         & veltry, velvec(kdir, natom)
      else
         print strfmt, natom, kdir, &
         & posold(kdir, natom), posvec(kdir, natom), posdif, &
         & veltry, velvec(kdir, natom)
      endif
   enddo
   enddo

   if ( int_stat == 0 ) then
      write( unit=funit, fmt=* ) ''
   else
      print *, ''
   endif
   close( unit=funit )

   posold = posvec

end subroutine Check_posvel
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
