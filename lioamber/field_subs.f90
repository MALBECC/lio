! This module contains routines for Field setups and calculations.
! Field shapes are as follows:
! 0 - Constant, 1 - Gaussian, 2 - Half Gaussian Up, 3 - Half Gaussian Down
module field_subs
   use shape_data, only: shape_iso, shape_aniso
   implicit none
   integer :: nfields_iso = 0, nfields_aniso = 0

   type, extends(shape_iso) :: field_iso
      real*8 :: amplitude = 0.0D0
   end type field_iso

   type, extends(shape_aniso) :: field_aniso
      real*8 :: amplitude(3) = 0.0D0
   end type field_aniso

   type(field_aniso), allocatable :: fields_aniso(:)
   type(field_iso)  , allocatable :: fields_iso(:)

   interface field_setup_iso
      module procedure field_setup_iso_full
      module procedure field_setup_iso_easy
   end interface field_setup_iso

contains

   ! Specifies field shape and parameters. the "easy" setup asumes default
   ! relations for some of the shape parameters. Both subroutines are mutually
   ! exclusive. In easy, the decay is applied arbtirarily.
   subroutine field_setup_iso_full(nfields_i, fld_shape, periodic, time_decay, &
                                   time_start, time_end, center, coord, amplit)
      implicit none
      integer, intent(in) :: nfields_i, fld_shape(nfields_i)
      logical, intent(in) :: periodic(nfields_i)
      real*8 , intent(in) :: time_decay(nfields_i), time_start(nfields_i), &
                             time_end(nfields_i), center(nfields_i),       &
                             coord(nfields_i,3), amplit
      integer :: icount

      nfields_iso = nfields_i
      if (allocated(fields_iso)) deallocate(fields_iso)
      allocate (fields_iso(nfields_iso))

      do icount = 1, nfields_iso
         fields_iso(icount)%shape_type = fld_shape(icount)
         fields_iso(icount)%periodic   = periodic(icount)
         fields_iso(icount)%time_decay = time_decay(icount)
         fields_iso(icount)%time_start = time_start(icount)
         fields_iso(icount)%time_end   = time_end(icount)
         fields_iso(icount)%center     = center(icount)
         fields_iso(icount)%coord(1:3) = coord(icount,1:3)
         fields_iso(icount)%amplitude  = amplit
      enddo
      return
   end subroutine field_setup_iso_full

   subroutine field_setup_iso_easy(nfields_i, fld_shape, time_end, coord)
      implicit none
      integer, intent(in) :: nfields_i, fld_shape(nfields_i)
      real*8 , intent(in) :: time_end(nfields_i), coord(nfields_i,3)
      integer :: icount

      nfields_iso = nfields_i
      if (allocated(fields_iso)) deallocate(fields_iso)
      allocate (fields_iso(nfields_iso))

      do icount = 1, nfields_iso
         fields_iso(icount)%shape_type = fld_shape(icount)
         fields_iso(icount)%periodic   = .false.
         fields_iso(icount)%time_start = 0.0D0
         fields_iso(icount)%time_end   = time_end(icount)
         fields_iso(icount)%coord(1:3) = coord(icount,1:3)

         select case (fld_shape(icount))
         case (0)
            fields_iso(icount)%time_decay = 0.0D0
            fields_iso(icount)%center     = 0.0D0
         case (1)
            fields_iso(icount)%time_decay = time_end(icount) /50.0D0
            fields_iso(icount)%center     = time_end(icount) / 2.0D0
         case (2)
            fields_iso(icount)%time_decay = time_end(icount) /25.0D0
            fields_iso(icount)%center     = time_end(icount)
         case (3)
            fields_iso(icount)%time_decay = time_end(icount) /25.0D0
            fields_iso(icount)%center     = 0.0D0
         case default
            write(*,*) "ERROR - Field: Wrong shape type for field ", icount, "."
            stop
         end select
      enddo
      return
   end subroutine field_setup_iso_easy

   subroutine field_calc()

   end subroutine field_calc

   subroutine field_finalize()
      implicit none
      if (allocated(fields_iso)) deallocate(fields_iso)
      if (allocated(fields_aniso)) deallocate(fields_aniso)
      return
   end subroutine field_finalize
end module field_subs
