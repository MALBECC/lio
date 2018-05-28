! This module contains shared variables for shape calculations.
! shape_type indicates the shape in time: 0 is constant, 1 is gaussian.
module shape_data
   private
   public :: shape_iso, shape_aniso

   type shape_iso
      integer :: shape_type = 0
      logical :: periodic   = .false.
      real*8  :: time_decay = 1.0D0
      real*8  :: time_start = 0.0D0
      real*8  :: time_end   = 0.0D0
      real*8  :: center     = 1.0D0
      real*8  :: coord(3)   = 0.0D0
      contains
         procedure, nopass :: calc_g
         procedure, nopass :: calc_c
         procedure, nopass :: calc_s
         procedure, pass   :: calculate => calc_iso
   end type shape_iso

   type shape_aniso
      integer :: shape_type(3) = 0
      logical :: periodic(3)   = .false.
      real*8  :: time_decay(3) = 1.0D0
      real*8  :: time_start(3) = 0.0D0
      real*8  :: time_end(3)   = 0.0D0
      real*8  :: center(3)     = 0.0D0
      real*8  :: coord(3)      = 0.0D0
      contains
         procedure, nopass :: calc_g
         procedure, nopass :: calc_c
         procedure, nopass :: calc_s
         procedure, pass   :: calculate => calc_aniso
   end type shape_aniso
contains

   ! Calculates the shape in one coordinate.
   ! Calculates a gaussian.
   function calc_g(in_time, f_coord, decay, center) result(shape_mag)
      implicit none
      real*8 , intent(in) :: in_time, f_coord, decay, center
      real*8  :: shape_mag, exponent, period

      shape_mag = 0.0D0
      exponent  = 0.0D0

      exponent  = - ( (in_time - center) / decay )**2
      shape_mag = f_coord * exp(exponent)

      return
   end function calc_g

   ! Calculates a sine.
   function calc_s(in_time, f_coord, period, phase) result(shape_mag)
      implicit none
      real*8, intent(in) :: in_time, f_coord, period, phase
      real*8 :: shape_mag, argument

      shape_mag = 0.0D0
      argument  = 0.0D0

      argument  = (in_time - phase) / period
      shape_mag = sin(argument)

      return
   end function calc_s

   ! Calculates a cosine.
   function calc_c(in_time, f_coord, period, phase) result(shape_mag)
      implicit none
      real*8, intent(in) :: in_time, f_coord, period, phase
      real*8 :: shape_mag, argument

      shape_mag = 0.0D0
      argument  = 0.0D0

      argument  = (in_time - phase) / period
      shape_mag = cos(argument)

      return
   end function calc_c

   subroutine calc_iso(this, time, shape_mag)
      implicit none
      class(shape_iso), intent(in)  :: this
      real*8          , intent(in)  :: time
      real*8          , intent(out) :: shape_mag(3)
      real*8  :: in_time, period
      integer :: icount, calc_coord_g

      if (.not.(this%periodic)) then
         if ((time.gt.this%time_end).or.(time.lt.this%time_start)) then
            shape_mag = 0.0D0
            return
         endif
         in_time = time
      else
         period  = this%time_end - this%time_start
         in_time = time  - floor(time/period) * period
      endif

      select case (this%shape_type)
         case (0)
            shape_mag = this%coord
         case (1, 2, 3)
            do icount = 1, 3
               shape_mag(icount) = this%calc_g(in_time, this%coord(icount), &
                                                this%time_decay, this%center)
            enddo
         case (4, 5)
            do icount = 1, 3
               shape_mag(icount) = this%calc_s(in_time, this%coord(icount), &
                                                this%time_decay, this%center)
            enddo
         case (6, 7)
            do icount = 1, 3
               shape_mag(icount) = this%calc_c(in_time, this%coord(icount), &
                                                this%time_decay, this%center)
            enddo
         case default
            write(*,*) "ERROR - shape. Wrong shape shape (calc_coord)."
      end select

      return
   end subroutine calc_iso

   subroutine calc_aniso(this, time, shape_mag)
      implicit none
      class(shape_aniso), intent(in)  :: this
      real*8            , intent(in)  :: time
      real*8            , intent(out) :: shape_mag(3)
      real*8  :: in_time(3), period
      integer :: icount, calc_coord_g

      do icount = 1, 3
         if((time.gt.this%time_end(icount)).and.(.not.(this%periodic(icount))))&
         then
            shape_mag(icount) = 0.0D0
         elseif (this%periodic(icount)) then
            period          = this%time_end(icount) - this%time_start(icount)
            in_time(icount) = time  - floor(time/period) * period
         else
            in_time(icount) = time
         endif

         if (time.le.this%time_end(icount)) then
            select case (this%shape_type(icount))
            case (0)
               shape_mag(icount) = this%coord(icount)
            case (1, 2, 3)
               shape_mag(icount) = this%calc_g(in_time(icount), &
                                   this%coord(icount), this%time_decay(icount),&
                                   this%center(icount))
            case (4, 5)
               shape_mag(icount) = this%calc_s(in_time(icount), &
                                   this%coord(icount), this%time_decay(icount),&
                                   this%center(icount))
            case (6,7)
               shape_mag(icount) = this%calc_c(in_time(icount), &
                                   this%coord(icount), this%time_decay(icount),&
                                   this%center(icount))
            case default
               write(*,*) "ERROR - shape. Wrong shape shape (calc_coord)."
            end select
         endif
      enddo

      return
   end subroutine calc_aniso
end module shape_data
