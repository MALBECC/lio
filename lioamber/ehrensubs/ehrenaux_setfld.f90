!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrenaux_setfld( current_time, elec_field )
!------------------------------------------------------------------------------!
!
! DESCRIPTION
! Time is in ps?fs?
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use ehrendata, &
   &only: eefld_ampx, eefld_ampy, eefld_ampz, eefld_wavelen                    &
       &, eefld_timegih, eefld_timegfh, eefld_timepos, eefld_timeamp

   implicit none
   real*8,intent(in)  :: current_time
   real*8,intent(out) :: elec_field(3)
   real*8  :: field_shape
   real*8  :: time_fact, time_dist, laser_freq
   logical :: apply_gaussian_gih, apply_gaussian_gfh

   real*8, parameter :: LIGHTSPEED_nm_fs = 299.792458d0
   real*8, parameter :: PI_x2 = 6.2831853071795864769d0
!
   field_shape = 1.0d0
!
!
!  GAUSSIAN SHAPE
   apply_gaussian_gih = (eefld_timegih) .and. (current_time < eefld_timepos)
   apply_gaussian_gfh = (eefld_timegfh) .and. (current_time > eefld_timepos)

   if (apply_gaussian_gih .or. apply_gaussian_gfh) then
      time_fact = (-1.0d0) / (eefld_timeamp)
      time_dist = (current_time - eefld_timepos)**2
      field_shape = field_shape * exp( (time_fact) * (time_dist) )
   endif
!
!
!  LASER SHAPE
   if (eefld_wavelen > 0.0d0) then
!      [ fs-1 ]                                         [long in nm]
      laser_freq = ( PI_x2 ) * ( LIGHTSPEED_nm_fs ) / ( eefld_wavelen )
      field_shape = field_shape * sin( laser_freq * current_time )
   endif
!
!
!  APPLY SHAPE
   elec_field(1) = eefld_ampx * field_shape
   elec_field(2) = eefld_ampy * field_shape
   elec_field(3) = eefld_ampz * field_shape
!
end subroutine ehrenaux_setfld
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
