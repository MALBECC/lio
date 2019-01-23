!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Time Dependent Analyze Program                                               !
! New TD analyze, based on U. Morzan and N. Foglia versions.                   !
! This program returns the oscillator strength's absolute value (in a.u.) vs   !
! energy (spectrum_eV) and wavelength (spectrum_nm).                           !
!                                                                              !
! Usage:                                                                       !
!                                                                              !
!                                                                              !
!                                                                              !
! December 2017 - F. Pedron                                                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

program tdanalyze

   implicit none
   ! External input variables
   ! Mu contains the time-dependent dipole moment (debyes) and its derivative
   ! dMu/dt afterwards. time_step contains the time step in fs. Damp is the
   ! spectrum damping factor (1/a.u. time units). field_x/y/z are the field
   ! strengths of the perturbation in each coordinate. n_steps is the total
   ! number of steps and lmin/lmax contain the spectrum limits.
   real*8, allocatable :: mu(:)
   real*8              :: time_step, damp, field
   integer             :: lmin, lmax, n_steps, column
   character(20)       :: input
   character(4)        :: lmax_i, lmin_i, damp_i, column_i

   ! Internal and auxiliar variables
   ! fti/ftr are the imaginary and real parts of the Fourier Transform of du/dt,
   ! nu is the frequency in 1/fs, time is the time in fs, and ene is the energy
   ! in eV.
   real*8, allocatable :: lambda(:), absorptivity(:)

   ! Reads input command line
   call logo()
   call get_command_argument(1, input)
   call get_command_argument(2, lmin_i)
   call get_command_argument(3, lmax_i)
   call get_command_argument(4, damp_i)
   call get_command_argument(5, column_i)

   ! Verifies and reads input file, allocating variables.
   call verify_input(input, lmin_i, lmax_i, damp_i, column_i, lmin, lmax, damp,&
                     column)
   call obtain_ns(n_steps, input)

   allocate (mu(n_steps), lambda(10*(lmax-lmin)), absorptivity(10*(lmax-lmin)))
   call read_input_file(input, n_steps, column, mu, time_step, field)

   ! Performs calculations and writes output.
   call calculate_spectrum(mu, field, time_step, n_steps, lmin, lmax, lambda,  &
                           absorptivity)
   call write_output(lmin, lmax, lambda, absorptivity)

contains

!%% VERIFY_INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine verify_input(input, lmin, lmax, damp, column, lmin_out, lmax_out,   &
                        damp_out, column_out)
   implicit none
   character(20), intent(in)    :: input
   character(4) , intent(inout) :: lmin, lmax, damp, column
   real*8       , intent(out)   :: damp_out
   integer      , intent(out)   :: column_out, lmin_out, lmax_out
   logical                      :: file_exists

   if (len_trim(input) .eq. 0) then
      write(*,*) "TDANALYZE needs a dipole moment input, please use:"
      write(*,*) "./tdanalyze3 dipole_moment_file"
      write(*,*) "or"
      write(*,*) "./tdanalyze3 dipole_moment_file lambda_min lambda_max &
                  damping_factor column_to_read"
      write(*,*)
      stop
   end if

   ! Verifies input file existence.
   inquire(file = input, exist = file_exists)
   if (.not.file_exists) then
      write(*,*) "File ", input, "does not exist, please verify file name."
      stop
   endif

   ! Checks the different input parameters.
   if (.not.(len_trim(lmin).gt.0)) then
      write(*,*) "Lambda minimum missing, please provide it (nm):"
      read(*,*) lmin
   endif

   if (.not.(len_trim(lmax).gt.0)) then
      write(*,*) "Lambda maximum missing, please provide it (nm):"
      read(*,*) lmax
   endif

   if (.not.(len_trim(damp).gt.0)) then
      write(*,*) "Damping factor missing, please provide it:"
      read(*,*) damp
   endif

   if (.not.((len_trim(column).gt.0).and.(len_trim(column).le.3))) then
      write(*,*) "Column to read outside ranks 1-3, please provide it:"
      read(*,*) column
   endif

   read(lmin  ,*) lmin_out
   read(lmax  ,*) lmax_out
   read(damp  ,*) damp_out
   read(column,*) column_out
   return
end subroutine verify_input
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% OBTAIN_NS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine obtain_ns(ns, input)
   implicit none
   integer  , intent(inout)  :: ns
   character(20), intent(in) :: input
   integer        :: ios
   real*8         :: any_line
   character(100) :: first_line

   ns = 0
   open (100, file=input)
   read (100,*) first_line
   do
      read (100,*,iostat=ios) any_line
      if (ios/=0) exit
      ns = ns + 1
   enddo
   close(100)

   write(*,201)  ns
201 format (1x,'Number of electron dynamic steps = ', I10)
endsubroutine obtain_ns
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% READ_INPUT_FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_input_file(input, ns, column, mu, time_step, field)
   implicit none
   character(20), intent(in)  :: input
   integer      , intent(in)  :: ns, column
   real*8       , intent(out) :: mu(ns), time_step, field
   real*8    :: t, mu_read(ns,3), field_read(3)
   integer   :: icount
   character :: hash(1)

   open (unit = 100, file = input)

   ! Reads initial hash [#], time step [fs], and field strengh of perturbation
   ! in the x, y and z axis [au]
   read(100,*) hash, time_step, field_read(1), field_read(2), field_read(3)
   field = field_read(column)

   ! Reads the dipole moment.
   do icount = 1, ns
     read(100,*) t, mu_read(icount,1), mu_read(icount,2), mu_read(icount,3)
     mu(icount) = mu_read(icount, column)
   enddo
   close(100)
end subroutine read_input_file
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% CALCULATE_SPECTRUM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine calculate_spectrum(mu, field_i, time_step, ns, lmin, lmax, lambda,  &
                              absorptivity)
   implicit none
   integer, intent(in)  :: ns, lmin, lmax
   real*8 , intent(in)  :: mu(ns), field_i, time_step
   real*8 , intent(out) :: lambda(10*(lmax-lmin)), absorptivity(10*(lmax-lmin))
   real*8  :: factor, field, dmu(ns-1), ftr, fti, ft_abs, nu, t
   integer :: icount, jcount

   ! Useful constants
   ! Light speed is in nm/fs
   double precision, parameter :: c  = 299.792458d0
   double precision, parameter :: pi = 3.1415926535897932384626433832795d0

   ! Unit corrections. Field is computed in V/nm.
   factor = 0.5292d0 / ( c * field_i * 2.542d0 * 0.02419d0 * 10d0 )
   field  = field_i * 514.220652d0

   ! Takes differences and makes dmu/dt
   do icount = 1, ns-1
      dmu(icount) = mu(icount+1) - mu(icount)
      dmu(icount) = dmu(icount)  / time_step
   enddo

   do icount = 1, 10*(lmax - lmin)
      lambda(icount) = dble((0.1d0*icount) + lmin)
      if ((mod(int(lambda(icount))-1,50).eq.0).and.(mod(icount-1,10) .eq. 0))  &
      then
         write(*,103) int(lambda(icount)), int(lambda(icount))+49
      end if

      ! Frequency in 1/fs
      nu     = c / lambda(icount)
      ft_abs = 0.0d0
      fti    = 0.0d0
      ftr    = 0.0d0
      t      = 0.0d0

      do jcount = 1, ns - 1
         t = t + time_step
         ! The real part of the transformation with a damping factor
         ftr = ftr + cos(2d0*pi*t*nu) * dmu(jcount) * exp(-t*damp)
         ! The imaginary part of the transformation with a damping factor
         fti = fti + sin(2d0*pi*t*nu) * dmu(jcount) * exp(-t*damp)
      enddo


      ! Takes absolute value and corrects units.
      ft_abs = ABS( DCMPLX(ftr, fti) )
!      write(*,*) ftr, fti, ft_abs
      ft_abs = ft_abs*time_step*4*pi
      absorptivity(icount) = ft_abs*factor
   enddo

   write (*,104) c*2.d0*time_step, 1.241462463D3/(c*2.d0*time_step)
   return
103 format (1x,"computing lambda from ",I5,1x,"to "I5)
104 format (1x, "Time step peak: ", F8.2, 1x, "nm, ", F10.2, 1x, "eV")
end subroutine calculate_spectrum
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% WRITE_OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_output(lmin, lmax, lambda, absorptivity)
   implicit none
   integer, intent(in) :: lmin, lmax
   real*8 , intent(in) :: lambda(10*(lmax-lmin)), absorptivity(10*(lmax-lmin))
   real*8  :: energy
   integer :: icount

   open (unit=100, file='espectro_eV')
   open (unit=101, file='espectro_nm')
   do icount = 1, 10*(lmax - lmin)
      energy = 1.241462463D3/lambda(icount)
      write (100, 100) energy        , absorptivity(icount)
      write (101, 100) lambda(icount), absorptivity(icount)
   enddo
   close(100)
   close(101)

   write (*,*)
   write (*,*) "TD Analyze finished."

   return
100 format (1x, E14.6, 2x, E14.6)
end subroutine write_output
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% LOGO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine logo()
   write(*,*)
   write(*,1200)
   write(*,1201)
   write(*,1202)
   write(*,1203)
   write(*,1204)
   write(*,1205)
   write(*,*)

1200 FORMAT(4x,"████████╗██████╗      █████╗ ███╗   ██╗ █████&
     ╗ ██╗     ██╗   ██╗███████╗███████╗")
1201 FORMAT(4x,"╚══██╔══╝██╔══██╗    ██╔══██╗████╗  ██║██&
     ╔══██╗██║     ╚██╗ ██╔╝╚══███╔╝██╔════╝")
1202 FORMAT(4x,"   ██║   ██║  ██║    ███████║██╔██╗ ██║███████&
     ║██║      ╚████╔╝   ███╔╝ █████╗  ")
1203 FORMAT(4x,"   ██║   ██║  ██║    ██╔══██║██║╚██╗██║██╔══██&
     ║██║       ╚██╔╝   ███╔╝  ██╔══╝  ")
1204 FORMAT(4x,"   ██║   ██████╔╝    ██║  ██║██║ ╚████║██║  ██║&
     ███████╗   ██║   ███████╗███████╗")
1205 FORMAT(4x,"   ╚═╝   ╚═════╝     ╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝╚═&
     ═════╝   ╚═╝   ╚══════╝╚══════╝")
end subroutine logo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

end program tdanalyze
