!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Time Dependent Analyze Program                                               !
! New TD analyze, based on U. Morzan and N. Foglia versions.                   !
! This program returns the oscillator strength's absolute value (in a.u.) vs   !
! energy (spectrum_eV) and wavelength (spectrum_nm).                           !
!                                                                              !
! Usage:                                                                       !
!   ./tdanalyze dipole_file lambda_min lamda_max damp_factor x_y_or_z          !
!                                                                              !
!                                                                              !
! March 2019 - F. Pedron                                                       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

program tdanalyze

   implicit none
   ! External input variables
   ! Mu contains the time-dependent dipole moment (debyes) and its derivative
   ! dMu/dt afterwards. time_step contains the time step in fs. Damp is the
   ! spectrum damping factor (1/a.u. time units). field_x/y/z are the field
   ! strengths of the perturbation in each coordinate. n_steps is the total
   ! number of steps and lmin/lmax contain the spectrum limits.
   real(kind=8), allocatable :: mu(:)
   character(len=5)  :: lmax_i, lmin_i, damp_i, column_i
   character(len=50) :: input
   integer           :: lmin, lmax, n_steps, column, err_stat
   real(kind=8)      :: time_step, damp, field

   ! Internal and auxiliar variables
   ! fti/ftr are the imaginary and real parts of the Fourier Transform of du/dt,
   ! nu is the frequency in 1/fs, time is the time in fs, and ene is the energy
   ! in eV.
   real(kind=8), allocatable :: lambda(:), absorptivity(:)

   ! Reads input command line
   call logo()
   call get_command_argument(1, input)
   call get_command_argument(2, lmin_i)
   call get_command_argument(3, lmax_i)
   call get_command_argument(4, damp_i)
   call get_command_argument(5, column_i)

   ! Verifies and reads input file, allocating variables.
   err_stat = 0
   call verify_input(input, lmin_i, lmax_i, damp_i, column_i, lmin, lmax, damp,&
                     column, err_stat)
   if (err_stat > 0) stop

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
                        damp_out, column_out, stat)
   implicit none
   character(len=50), intent(in)    :: input
   character(len=5) , intent(inout) :: lmin, lmax, damp, column
   real(kind=8)     , intent(out)   :: damp_out
   integer          , intent(out)   :: column_out, lmin_out, lmax_out, stat

   logical :: file_exists

   if (len_trim(input) .eq. 0) then
      write(*,'(A)') "TDANALYZE needs a dipole moment input, please use:"
      write(*,'(A)') "  ./tdanalyze dipole_moment_file"
      write(*,'(A)') "or"
      write(*,'(A)') "  ./tdanalyze dipole_moment_file lambda_min lambda_max &
                     & damping_factor column_to_read"
      write(*,*)
      stat = 1
      return
   end if

   ! Verifies input file existence.
   inquire(file = input, exist = file_exists)
   if (.not. file_exists) then
      write(*,'(A)') "File ", trim(input), &
                     "does not exist, please verify file name."
      stat = 1
      return
   endif

   ! Checks the different input parameters.
   if (.not. (len_trim(lmin) > 0)) then
      write(*,'(A)') "Lambda minimum missing, please provide it (nm):"
      read(*,*) lmin
   endif

   if (.not. (len_trim(lmax) > 0)) then
      write(*,'(A)') "Lambda maximum missing, please provide it (nm):"
      read(*,*) lmax
   endif

   if (.not. (len_trim(damp) > 0)) then
      write(*,'(A)') "Damping factor missing, please provide it:"
      read(*,*) damp
   endif

   if (.not. ((len_trim(column) > 0) .and. (len_trim(column) <= 3))) then
      write(*,'(A)') "Column to read outside ranks 1-3, please provide it:"
      read(*,*) column
   endif

   read(lmin  ,*) lmin_out
   read(lmax  ,*) lmax_out
   read(damp  ,*) damp_out
   read(column,*) column_out
end subroutine verify_input
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% OBTAIN_NS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine obtain_ns(ns, input)
   implicit none
   character(len=50), intent(in)    :: input
   integer          , intent(inout) :: ns

   character(len=100) :: first_line
   character(len=14)  :: fmt
   integer            :: ios
   real(kind=8)       :: any_line

   ns = 0
   open(100, file = input)
   read(100,*) first_line
   do
      read(100,*,iostat = ios) any_line
      if (ios /= 0) exit
      ns = ns + 1
   enddo
   close(100)

   fmt = '(1x, A35, I10)'
   write(*,fmt) 'Number of electron dynamic steps = ', ns
endsubroutine obtain_ns
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% READ_INPUT_FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_input_file(input, ns, column, mu, time_step, field)
   implicit none
   character(len=20), intent(in)  :: input
   integer          , intent(in)  :: ns, column
   real(kind=8)     , intent(out) :: mu(ns), time_step, field

   character(len=1) :: hash
   integer          :: icount
   real(kind=8)     :: t, mu_read(ns,3), field_read(3)

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
   integer     , intent(in)  :: ns, lmin, lmax
   real(kind=8), intent(in)  :: mu(ns), field_i, time_step
   real(kind=8), intent(out) :: lambda(10*(lmax-lmin)), &
                                absorptivity(10*(lmax-lmin))
   
   character(len=38) :: fmt
   integer           :: icount, jcount
   real(kind=8)      :: factor, dmu(ns-1), ftr, fti, ft_abs, nu, t

   ! Useful constants
   ! Light speed is in nm/fs
   real(kind=8), parameter :: C  = 299.792458d0
   real(kind=8), parameter :: PI = 3.1415926535897932384626433832795d0

   ! Unit corrections. Field is computed in V/nm.
   factor = 0.5292d0 / ( C * field_i * 2.542d0 * 0.02419d0 * 10d0 )

   ! Takes differences and makes dmu/dt
   do icount = 1, ns-1
      dmu(icount) = mu(icount+1) - mu(icount)
      dmu(icount) = dmu(icount) / time_step
   enddo

   fmt = '(1x,A22,I5,1x,A3,I5)'
   do icount = 1, 10*(lmax - lmin)
      lambda(icount) = dble((0.1d0*icount) + lmin)
      if ((mod(int(lambda(icount))-1,50) == 0) .and. (mod(icount-1,10) == 0)) &
      then
         write(*,fmt) "computing lambda from ", int(lambda(icount)), "to ", &
                      int(lambda(icount))+49
      endif

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
      ft_abs = ft_abs * time_step * 4d0 * pi
      absorptivity(icount) = ft_abs * factor
   enddo

   fmt = '(1x, A16, F8.2, 1x, A4, F10.2, 1x, A2)'
   write (*,fmt) "Time step peak: ", C * 2.d0 * time_step, "nm, ", &
                 1.241462463D3 / (C * 2.d0 * time_step), "eV"
end subroutine calculate_spectrum
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% WRITE_OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_output(lmin, lmax, lambda, absorptivity)
   implicit none
   integer     , intent(in) :: lmin, lmax
   real(kind=8), intent(in) :: lambda(10*(lmax-lmin)), &
                               absorptivity(10*(lmax-lmin))

   character(len=22) :: fmt
   integer           :: icount
   real(kind=8)      :: energy

   fmt = '(1x, E14.6, 2x, E14.6)'
   open (unit=100, file='espectro_eV')
   open (unit=101, file='espectro_nm')
   do icount = 1, 10*(lmax - lmin)
      energy = 1.241462463D3/lambda(icount)
      write (100, fmt) energy        , absorptivity(icount)
      write (101, fmt) lambda(icount), absorptivity(icount)
   enddo
   close(100)
   close(101)

   write (*,*)
   write (*,'(A)') "TD Analyze finished."

end subroutine write_output
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% LOGO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine logo()
   character(len=7) :: fmt = '(4x, A)'

   write(*,*)
   write(*,fmt) "████████╗██████╗      █████╗ ███╗   ██╗ █████&
   &╗ ██╗   ██╗   ██╗███████╗███████╗"
   write(*,fmt) "╚══██╔══╝██╔══██╗    ██╔══██╗████╗  ██║██&
   &╔══██╗██║   ╚██╗ ██╔╝╚══███╔╝██╔════╝" 
   write(*,fmt) "   ██║   ██║  ██║    ███████║██╔██╗ ██║███████&
   &║██║    ╚████╔╝   ███╔╝ █████╗  " 
   write(*,fmt) "   ██║   ██║  ██║    ██╔══██║██║╚██╗██║██╔══██&
   &║██║     ╚██╔╝   ███╔╝  ██╔══╝  " 
   write(*,fmt) "   ██║   ██████╔╝    ██║  ██║██║ ╚████║██║  ██║&
   &███████╗ ██║  ███████╗███████╗" 
   write(*,fmt) "   ╚═╝   ╚═════╝     ╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝╚═&
   &═════╝ ╚═╝   ╚══════╝╚══════╝" 
   write(*,*)

end subroutine logo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

end program tdanalyze
