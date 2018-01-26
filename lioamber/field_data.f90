! This module contains shared variables for field calculations using intfld.
! Shape indicates the field shape in time: 0 is constant, 1 is gaussian,  2 is
! half gaussian up, 3 is half gaussian down.
module field_data
   private
   public :: field, field_asym
   public :: nfields, nfields_asym, all_fields, all_asym_fields

   type field
      integer :: fld_shape  = 0
      logical :: periodic   = .false.
      real*8  :: time_decay = 1.0D0
      real*8  :: time_start = 0.0D0
      real*8  :: time_end   = 0.0D0
      real*8  :: center     = 1.0D0
      real*8  :: amplitude  = 1.0D0
      real*8  :: coord(3)   = 0.0D0
      contains
         procedure, nopass :: calc_g
         procedure, pass   :: calculate => calc_sym
   end type field

   type field_asym
      integer :: fld_shape(3)  = 0
      logical :: periodic(3)   = .false.
      real*8  :: time_decay(3) = 1.0D0
      real*8  :: time_start(3) = 0.0D0
      real*8  :: time_end(3)   = 0.0D0
      real*8  :: center(3)     = 0.0D0
      real*8  :: amplitude(3)  = 1.0D0
      real*8  :: coord(3)      = 0.0D0
      contains
         procedure, nopass :: calc_g
         procedure, pass   :: calculate => calc_asym
   end type field_asym

   integer :: nfields = 0
   integer :: nfields_asym = 0
   type(field)     , allocatable :: all_fields(:)
   type(field_asym), allocatable :: all_asym_fields(:)
contains

   ! Calculates the field in one coordinate.
   function calc_g(in_time, f_coord, decay, center) result(field_mag)
      implicit none
      real*8 , intent(in) :: in_time, f_coord, decay, center
      real*8  :: field_mag, exponent, period

      field_mag = 0.0D0
      exponent  = 0.0D0

      exponent  = - ( (in_time - center) / decay )**2
      field_mag = f_coord * exp(exponent)

      return
   end function calc_g

   subroutine calc_sym(this, time, field_mag)
      implicit none
      class(field), intent(in)  :: this
      real*8     , intent(in)  :: time
      real*8     , intent(out) :: field_mag(3)
      real*8  :: in_time, period
      integer :: icount, calc_coord_g

      if ( (time.gt.this%time_end).and.(.not.(this%periodic)) ) then
         field_mag = 0.0D0
         return
      elseif (this%periodic) then
         period  = this%time_end - this%time_start
         in_time = time  - floor(time/period) * period
      else
         in_time = time
      endif

      select case (this%fld_shape)
         case (0)
            field_mag = this%coord
         case (1, 2, 3)
            do icount = 1, 3
               field_mag(icount) = this%calc_g(in_time, this%coord(icount), &
                                                this%time_decay, this%center)
            enddo
         case default
            write(*,*) "ERROR - Field. Wrong field shape (calc_coord)."
      end select

      return
   end subroutine calc_sym

   subroutine calc_asym(this, time, field_mag)
      implicit none
      class(field_asym), intent(in)  :: this
      real*8          , intent(in)  :: time
      real*8          , intent(out) :: field_mag(3)
      real*8  :: in_time(3), period
      integer :: icount, calc_coord_g

      do icount = 1, 3
         if((time.gt.this%time_end(icount)).and.(.not.(this%periodic(icount))))&
         then
            field_mag(icount) = 0.0D0
         elseif (this%periodic(icount)) then
            period          = this%time_end(icount) - this%time_start(icount)
            in_time(icount) = time  - floor(time/period) * period
         else
            in_time(icount) = time
         endif

         select case (this%fld_shape(icount))
         case (0)
            if (time.le.this%time_end(icount)) then
               field_mag(icount) = this%coord(icount)
            endif
         case (1, 2, 3)
            if (time.le.this%time_end(icount)) then
               field_mag(icount) = this%calc_g(in_time(icount), &
                                   this%coord(icount), this%time_decay(icount),&
                                   this%center(icount))
            endif
         case default
            write(*,*) "ERROR - Field. Wrong field shape (calc_coord)."
         end select
      enddo

      return
   end subroutine calc_asym
end module field_data
