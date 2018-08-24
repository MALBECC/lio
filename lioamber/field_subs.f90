! FIELD_SUBS - FIELD_DATA
!
! This module contains routines for Field setups and calculations.
! Field shapes are as follows:
! 0 - Constant, 1 - Gaussian    , 2 - Half Gaussian Up, 3 - Half Gaussian Down
! 4 - Sine    , 5 - Shifted Sine, 6 - Cosine          , 7 - Shifted Cosine
! In both sine and cosine "decay" means "period" and "center" means "phase".
! Field ISO functions share the same f(t) in x, y and z coordinates (i.e
! x=A*f(t), y=B*g(t), z=C*h(t)), while ANISO functions are for example
! x=A*f(t), y=B*g(t), z=C*h(t).
!
! Input file reading when nfields > 0:
! For ISO fields, the input file should be an ASCII file with information
! pertaining each field in a separate line, containing (in this order) shape,
! periodic, time start, time end, decay/period, center/phase,  magnitude in x,
! magnitude in y, magnitude in z. For example:
!
! 1 f 0.0D0 0.2D0 4.0D-3 0.1D0 0.0D0 0.0D0 0.05
! This produces the standart TD gaussian perturbation in a Z axis.
!
! For ANISO fields, each line should contain the data for each coordinate of a
! field, meaning there are 3*n_fields number of lines (even if a default is
! desired, it should be manually set). Similar to the previous case, the line
! should contain shape, periodic, time start, time end, decay/period,
! center/phase and magnitude.

module field_data
   use shape_data, only: shape_iso, shape_aniso
   implicit none

   type, extends(shape_iso) :: field_iso
      real*8 :: field_g = 0.0D0
   end type field_iso
   type, extends(shape_aniso) :: field_aniso
      real*8 :: field_g = 0.0D0
   end type field_aniso

   logical :: field         = .false.
   integer :: nfields_iso   = 0
   integer :: nfields_aniso = 0
   real*8  :: chrg_sq       = 0.0D0
   real*8  :: epsilon       = 1.0D0
   real*8  :: a0            = 1000.0D0
   real*8  :: fx            = 0.0D0
   real*8  :: fy            = 0.0D0
   real*8  :: fz            = 0.0D0
   character*20 :: field_iso_file   = "field.in"
   character*20 :: field_aniso_file = "fieldaniso.in"
   type(field_aniso), allocatable :: fields_aniso(:)
   type(field_iso)  , allocatable :: fields_iso(:)
end module field_data

module field_subs
   use field_data, only: field_iso, field_aniso
   implicit none

   interface field_calc
      module procedure field_calc
   end interface field_calc

   interface read_fields
      module procedure read_fields
   end interface read_fields

   interface field_setup_iso
      module procedure field_setup_iso_full
      module procedure field_setup_iso_easy
   end interface field_setup_iso

   interface field_setup_aniso
      module procedure field_setup_aniso_full
      module procedure field_setup_aniso_easy
   end interface field_setup_aniso
contains

   ! Specifies field shape and parameters. the "easy" setup asumes default
   ! relations for some of the shape parameters. Both subroutines are mutually
   ! exclusive. In easy, the decay is applied arbtirarily.
   subroutine field_setup_iso_full(nfields_i, fld_shape, periodic, time_decay, &
                                   time_start, time_end, center, coord)
      use field_data, only: nfields_iso, fields_iso
      implicit none
      integer, intent(in) :: nfields_i, fld_shape(nfields_i)
      logical, intent(in) :: periodic(nfields_i)
      real*8 , intent(in) :: time_decay(nfields_i), time_start(nfields_i), &
                             time_end(nfields_i), center(nfields_i),       &
                             coord(nfields_i,3)
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
         fields_iso(icount)%field_g    = 1.0D0
      enddo
      return
   end subroutine field_setup_iso_full

   subroutine field_setup_iso_easy(nfields_i, fld_shape, time_end, coord)
      use field_data, only: nfields_iso, fields_iso
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
         fields_iso(icount)%field_g    = 1.0D0

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
         case (4)
            fields_iso(icount)%time_decay = time_end(icount)
            fields_iso(icount)%center     = 0.0D0
         case (5)
            fields_iso(icount)%time_decay = time_end(icount)
            fields_iso(icount)%center     = time_end(icount) / 4.0D0
         case (6)
            fields_iso(icount)%time_decay = time_end(icount)
            fields_iso(icount)%center     = 0.0D0
         case (7)
            fields_iso(icount)%time_decay = time_end(icount)
            fields_iso(icount)%center     = time_end(icount) / 4.0D0
         case default
            write(*,*) "ERROR - Field: Wrong shape type for field ", icount, "."
            stop
         end select
      enddo
      return
   end subroutine field_setup_iso_easy

   subroutine field_setup_aniso_full(nfields_i, fld_shape, period, time_decay, &
                                     time_start, time_end, center,  coord)
      use field_data, only: nfields_aniso, fields_aniso
      implicit none
      integer, intent(in) :: nfields_i, fld_shape(nfields_i,3)
      logical, intent(in) :: period(nfields_i,3)
      real*8 , intent(in) :: time_decay(nfields_i,3), time_start(nfields_i,3), &
                             time_end(nfields_i,3), center(nfields_i,3),       &
                             coord(nfields_i,3)
      integer :: icount

      nfields_aniso = nfields_i
      if (allocated(fields_aniso)) deallocate(fields_aniso)
      allocate (fields_aniso(nfields_aniso))

      do icount = 1, nfields_aniso
         fields_aniso(icount)%shape_type(1:3) = fld_shape(icount,1:3)
         fields_aniso(icount)%periodic(1:3)   = period(icount,1:3)
         fields_aniso(icount)%time_decay(1:3) = time_decay(icount,1:3)
         fields_aniso(icount)%time_start(1:3) = time_start(icount,1:3)
         fields_aniso(icount)%time_end(1:3)   = time_end(icount,1:3)
         fields_aniso(icount)%center(1:3)     = center(icount,1:3)
         fields_aniso(icount)%coord(1:3)      = coord(icount,1:3)
         fields_aniso(icount)%field_g         = 1.0D0
      enddo
      return
   end subroutine field_setup_aniso_full

   subroutine field_setup_aniso_easy(nfields_i, fld_shape, tm_end, coord)
      use field_data, only: nfields_aniso, fields_aniso
      implicit none
      integer, intent(in) :: nfields_i, fld_shape(nfields_i,3)
      real*8 , intent(in) :: tm_end(nfields_i,3), coord(nfields_i,3)
      integer :: icount, crd

      nfields_aniso = nfields_i
      if (allocated(fields_aniso)) deallocate(fields_aniso)
      allocate (fields_aniso(nfields_aniso))

      do icount = 1, nfields_aniso
         fields_aniso(icount)%shape_type(1:3) = fld_shape(icount,1:3)
         fields_aniso(icount)%periodic(1:3)   = .false.
         fields_aniso(icount)%time_start(1:3) = 0.0D0
         fields_aniso(icount)%time_end(1:3)   = tm_end(icount,1:3)
         fields_aniso(icount)%coord(1:3)      = coord(icount,1:3)
         fields_aniso(icount)%field_g         = 1.0D0

         do crd = 1, 3
            select case (fld_shape(icount, crd))
            case (0)
               fields_aniso(icount)%time_decay(crd) = 0.0D0
               fields_aniso(icount)%center(crd)     = 0.0D0
            case (1)
               fields_aniso(icount)%time_decay(crd) = tm_end(icount,crd) /50.0D0
               fields_aniso(icount)%center(crd)     = tm_end(icount,crd) / 2.0D0
            case (2)
               fields_aniso(icount)%time_decay(crd) = tm_end(icount,crd) /25.0D0
               fields_aniso(icount)%center(crd)     = tm_end(icount,crd)
            case (3)
               fields_aniso(icount)%time_decay(crd) = tm_end(icount,crd) /25.0D0
               fields_aniso(icount)%center(crd)     = 0.0D0
            case (4)
               fields_aniso(icount)%time_decay(crd) = tm_end(icount,crd)
               fields_aniso(icount)%center(crd)     = 0.0D0
            case (5)
               fields_aniso(icount)%time_decay(crd) = tm_end(icount,crd)
               fields_aniso(icount)%center(crd)     = tm_end(icount,crd) / 4.0D0
            case (6)
               fields_aniso(icount)%time_decay(crd) = tm_end(icount,crd)
               fields_aniso(icount)%center(crd)     = 0.0D0
            case (7)
               fields_aniso(icount)%time_decay(crd) = tm_end(icount,crd)
               fields_aniso(icount)%center(crd)     = tm_end(icount,crd) / 4.0D0

            case default
               write(*,*) "ERROR - Field: Wrong shape type for field ", icount,&
                          "."
               stop
            end select
         enddo
      enddo

      return
   end subroutine field_setup_aniso_easy

   ! Intended to use as default or retrocompatibility with old way of
   ! setting fields.
   subroutine field_setup_old(pert_time, fld_shape, fld_x, fld_y, fld_z)
      use field_data, only: nfields_iso, nfields_aniso
      implicit none
      integer, intent(in) :: fld_shape
      real*8 , intent(in) :: fld_x, fld_y, fld_z, pert_time
      real*8  :: coord(3) , time_end(1)
      integer :: nfields_i, shape_i(1)

      nfields_i   = 1
      shape_i(1)  = fld_shape
      time_end(1) = pert_time
      coord(1)    = fld_x
      coord(2)    = fld_y
      coord(3)    = fld_z
      call field_setup_iso_easy(nfields_i, shape_i, time_end, coord)

      return
   end subroutine field_setup_old

   ! Performs field calculation adding all field magnitudes to x, y and z
   ! coordinates.
   subroutine field_calc_all(fx, fy, fz, time)
      use field_data, only: nfields_iso, nfields_aniso, fields_iso, fields_aniso
      real*8, intent(in)  :: time
      real*8, intent(out) :: fx, fy, fz
      real*8              :: fld_temp(3)
      type(field_iso)     :: f_temp
      integer             :: icount

      fx = 0.0D0; fy = 0.0D0; fz = 0.0D0
      if (allocated(fields_iso)) then
         do icount = 1, nfields_iso
            fld_temp = 0.0D0
            call fields_iso(icount)%calculate(time, fld_temp)
            fx = fx + fld_temp(1)
            fy = fy + fld_temp(2)
            fz = fz + fld_temp(3)
         enddo
      endif

      if (allocated(fields_aniso)) then
         do icount = 1, nfields_aniso
            fld_temp = 0.0D0
            call fields_aniso(icount)%calculate(time, fld_temp)
            fx = fx + fld_temp(1)
            fy = fy + fld_temp(2)
            fz = fz + fld_temp(3)
         enddo
      endif

      return
   end subroutine field_calc_all

   subroutine field_calc(energ, time, Fmat, Fmat_b, r, d, Iz, natom, ntatom, &
                         opshell)
      use faint_cpu , only: intfld
      use field_data, only: chrg_sq, epsilon, a0
      implicit none
      integer, intent(in)             :: natom, ntatom, Iz(ntatom)
      logical, intent(in)             :: opshell
      double precision, intent(in)    :: time, r(ntatom,3), d(natom,natom)
      double precision, intent(inout) :: energ, Fmat(:), Fmat_b(:)
      double precision :: dipxyz(3), g, factor, Fx, Fy, Fz, tol

      Fx     = 0.0D0  ; Fy     = 0.0D0
      Fz     = 0.0D0  ; energ  = 0.0D0
      dipxyz = 0.0D0  ; g      = 1.00D0
      factor = 2.54D0 ; tol    = 1.00D-16

      call dip(dipxyz)
      call field_calc_all(Fx, Fy, Fz, time)
      if ((abs(Fx).lt.tol) .and. (abs(Fy).lt.tol) .and. (abs(Fz).lt.tol)) return
      call intfld(Fmat, Fmat_B, r, d, Iz, natom, ntatom, opshell, g, Fx, Fy, Fz)
      energ = - g * (Fx*dipxyz(1) + Fy*dipxyz(2) + Fz*dipxyz(3)) / factor -   &
                0.50D0 * (1.0D0 - 1.0D0/epsilon) * chrg_sq/a0
      return
   end subroutine field_calc

   ! Reads field from input files.
   subroutine read_iso_fields()
      use field_data, only: nfields_iso, field_iso_file
      implicit none
      integer, allocatable :: fld_shape(:)
      logical, allocatable :: periodic(:)
      real*8 , allocatable :: time_decay(:), time_start(:), time_end(:),  &
                              center(:), coord(:,:)
      logical :: file_exists
      integer :: icount, ios

      if (nfields_iso.lt.1) return
      allocate (fld_shape(nfields_iso) , periodic(nfields_iso)  , &
                time_start(nfields_iso), time_decay(nfields_iso), &
                time_end(nfields_iso)  , center(nfields_iso)    , &
                coord(nfields_iso, 3))

      inquire(file = field_iso_file, exist = file_exists)
      if (.not.file_exists) then
         stop "ERROR - field_iso file not found, but nfields_iso > 0."
      endif

      open(unit=314, file=field_iso_file, iostat=ios)
      do icount = 1, nfields_iso
         read(314,*) fld_shape(icount), periodic(icount)  , time_start(icount),&
                     time_end(icount) , time_decay(icount), center(icount)    ,&
                     coord(icount,1)  , coord(icount,2)   , coord(icount,3)
      enddo
      close(314)
      call field_setup_iso_full(nfields_iso, fld_shape, periodic, time_decay,  &
                                time_start , time_end , center  ,  coord)
      deallocate (fld_shape, periodic, time_start, time_decay, time_end, &
                  center, coord)
      return
   end subroutine read_iso_fields

   subroutine read_aniso_fields()
      use field_data, only: nfields_aniso, field_aniso_file
      implicit none
      integer, allocatable :: fld_shape(:,:)
      logical, allocatable :: periodic(:,:)
      real*8 , allocatable :: time_decay(:,:), time_start(:,:), time_end(:,:), &
                              center(:,:), coord(:,:)
      logical :: file_exists
      integer :: icount, jcount, ios

      if (nfields_aniso.lt.1) return
      allocate (fld_shape(nfields_aniso,3) , periodic(nfields_aniso,3)  , &
                time_start(nfields_aniso,3), time_decay(nfields_aniso,3), &
                time_end(nfields_aniso,3)  , center(nfields_aniso,3)    , &
                coord(nfields_aniso, 3))

      inquire(file = field_aniso_file, exist = file_exists)
      if (.not.file_exists) then
         stop "ERROR - field_iso file not found, but nfields_iso > 0."
      endif

      open(unit=314, file=field_aniso_file, iostat=ios)
      do icount = 1, nfields_aniso
      do jcount = 1, 3
         read(314,*) fld_shape(icount,jcount) , periodic(icount,jcount), &
                     time_start(icount,jcount), time_end(icount,jcount), &
                     time_decay(icount,jcount), center(icount,jcount)  , &
                     coord(icount,jcount)
      enddo
      enddo
      close(314)
      call field_setup_aniso_full(nfields_aniso, fld_shape, periodic, &
                                  time_decay, time_start , time_end , center, &
                                  coord)
      deallocate (fld_shape, periodic, time_start, time_decay, time_end, &
                  center, coord)
      return
   end subroutine read_aniso_fields

   subroutine read_fields()
      call read_iso_fields()
      call read_aniso_fields()
   end subroutine read_fields

   ! Deallocates field arrays.
   subroutine field_finalize()
      use field_data, only: fields_iso, fields_aniso
      implicit none
      if (allocated(fields_iso)) deallocate(fields_iso)
      if (allocated(fields_aniso)) deallocate(fields_aniso)
      return
   end subroutine field_finalize
end module field_subs
