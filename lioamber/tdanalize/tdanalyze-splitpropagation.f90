!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Time Dependent Analyze Program                                                                                                   !
!                                                                                                                                  !
! modification for tdanalize1 of U. Morzan made by Nick                                                                            !
! This program return absolute value of osc. stenght in a.u. vs energy (espectro_eV) and vs wavelength (espectro_nm)               !
!                                                                                                                                  !
! V 1.00 April  2016                                                                                                               !
!                                                                                                                                  !
! Nicolas Foglia                                                                                                                   !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!call this program as:                                                                                                             !
! ./tdanalyze3 dipole_moment_input                                                                                                 !
! or                                                                                                                               !
! ./tdanalyze3 dipole_moment_input lambda_min lambda_max damping_factor                                                            !
!                                                                                                                                  !
! dipole_moment_input:                                                                                                             !
!   ts  NCO  field                                                                                                                 !
!   trash-line                                                                                                                     !
!   time1   dipole1.1 dipole1.2 dipole1.3 ...  dipole1.NCO                                                                         !
!   time2   dipole2.1 dipole2.2 dipole2.3 ...  dipole2.NCO                                                                         !
!   ...                                                                                                                            !
!   ...                                                                                                                            !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


      PROGRAM TDANALYZE
!Variables description:
!mu at first have dipole moment [Debyes] after is transform into du/dt [Debyes/fs]
!fti,ftr imaginary and real part of fourier transform of du/dt
!ts = time step (fs)
!nu = frequency computed in 1/fs
!t = time[fs]
!ene = energy [eV]
!damp damping factor [(a.u. time)^-1]
!field strengh of perturbation [a.u.]
!ns = numbers of dinamic steps
!lmin,lmax lambda for spectrum calculation
!NCO = number of occupied orbitals, NCO=1 for not split propagation

       implicit none
       double precision, allocatable :: mu(:,:),fti(:),ftr(:)
       double precision :: ts,nu,t,ene,damp,field
       integer :: ns,lmin,lmax,NCO

       double precision :: factor !factor for correct osc. strenght to au
       double precision :: lambda !auxiliar
       character(20)  :: trash !auxiliar
       integer :: i,j,l,extra !auxiliar


       double precision, parameter :: c=299.792458d0 !light speed in nm/fs
       double precision, parameter :: pi=3.1415926535897932384626433832795d0

       character(20) input! dipole_moment_input
       character(len=24) :: formatout !variable format for split propagation output
       character(4) lminr,lmaxr,dampr !input variables for lambda min, max and damping factor
       LOGICAL :: existfile


       call get_command_argument(1,input)
       call get_command_argument(2,lminr)
       call get_command_argument(3,lmaxr)
       call get_command_argument(4,dampr)

       ns=0
       call logo()


       if (len_trim(input) .eq. 0) then
        write(*,*) "TD ANALYZE needs a dipole_moment_input, please use:"
        write(*,*) "./tdanalyze3 dipole_moment_input"
        write(*,*) "or"
        write(*,*) "./tdanalyze3 dipole_moment_input lambda_min &
        lambda_max damping_factor "
        write(*,*)
        STOP
       end if

       INQUIRE(FILE=INPUT, EXIST=existfile)
       IF ( .NOT. existfile) THEN !verifica la existencia del archivo de parametros
        WRITE(*,*) input, "does not exist, please verify file name"
        STOP
       END IF

       if (len_trim(lminr).gt.0 .and. len_trim(lmaxr).gt.0 .and. len_trim(dampr) .gt. 0) then
        read(lminr,*) lmin
        read(lmaxr,*) lmax
        read(dampr,*) damp
       else
        write(*,*)
        write(*,*)"please write:"
        write(*,*)'lambda min(nm),lambda max(nm), damp(a.u.)'
        read(*,*) lmin, lmax, damp
        write(*,*)
       endif

       call obtain_ns(ns,input) !obtain number of steps in input

       open (unit=1,file=input)
!       open (unit=2,file=output)
       open (unit=3,file='espectro_eV')
!       open (unit=4, file='sq_espectro')
       open (unit=5, file='espectro_nm')

       read(1,*) ts, NCO, field ! time step [fs],number of occupied orbitals,field strengh of perturbation [au]
       read(1,*) trash !not used


       if (NCO .eq. 1) then !special case for not split propagation
        allocate (mu(1:ns,NCO),fti(NCO),ftr(NCO))
        extra=0
        write (formatout, "('(1x,E14.6,',I3.3,'(2x,E14.6))')") NCO !define formato de salida
       else
        allocate (mu(1:ns,NCO+1),fti(NCO+2),ftr(NCO+2))
        extra=1
        write (formatout, "('(1x,E14.6,',I3.3,'(2x,E14.6))')") NCO+2 !define formato de salida
       end if

       factor=0.5292d0/(c*field*2.542d0*0.02419d0*10d0)
       field=field*514.220652d0 !field is computed in V/nm

       do i = 1, ns
        read(1,*) t, mu(i,:) !reading dipole moment
       enddo

       DO l=1,NCO+extra
        do i = 1, ns-1
         mu(i,l) = mu(i+1,l) - mu(i,l)  !takes differences
         mu(i,l)=mu(i,l)/ts           !mu now takes du/dt
        enddo
       END DO

       fti=0.d0
       ftr=0.d0

       do i = 1, 10*(lmax - lmin)
        lambda = dble((0.1d0*i) + lmin)
        if (mod(int(lambda)-1,50) .eq. 0 .and. mod(i-1,10) .eq. 0) then
         write(*,103) int(lambda),int(lambda)+49
        end if

        nu=c/lambda !frequency computed in 1/fs
        fti = 0.0d0
        ftr = 0.0d0
        t=0.0

        do j = 1,ns-1
         t = t + ts
         do l=1,NCO+extra
          ftr(l) = ftr(l) + cos(2d0*pi*t*nu) * mu(j,l)*exp(-t*damp) ! the real part of the transformation with a damping factor
          fti(l) = fti(l) + sin(2d0*pi*t*nu) * mu(j,l)*exp(-t*damp) ! the imaginary part of the transformation with a damping factor
         end do
        enddo


        ene = 1.241462463D3/lambda

        IF (extra .eq. 1) Then
         DO l=1,NCO+1
          ftr(NCO+2)=ftr(NCO+2)+ftr(l)
          fti(NCO+2)=fti(NCO+2)+fti(l)
         END DO
        end if


        DO l=1,NCO+2*extra
         ftr(l)=ABS(DCMPLX(ftr(l),fti(l)))  !using abs for take spectrun in correct fase
         ftr(l)=ftr(l)*ts*4*pi!a/(4*pi*pi)
        ENDDO

        write (3,formatout) ene, ftr*factor
        write (5,formatout) lambda, ftr*factor

!           write (2,100) lambda, ene, ftr, fti
       enddo

       write (*,*)
       write (*,104) c*2.d0*ts, 1.241462463D3/(c*2.d0*ts)
       write (*,*)
       write (*,*) "calculations END"
       close(1)
!      close(2)
       close(3)
!       close(4)
       close(5)

100       format (1x,E14.6,2x,E14.6,2x,E14.6,2x,E14.6)
101       format (1x,E14.6,2x,E14.6)
102       format (1x,E14.6,30(2x,E14.6))
103       format (1x,"computing lambda from ",I5,1x,"to "I5)
104       format (1x,"time step peak ", F8.2,1x,"nm, ", F10.2,1x,"eV")

       contains
       subroutine obtain_ns(ns,input)
          integer, intent(inout):: ns
          integer ::ios
          double precision ::  whatever
          character(LEN=100):: command
          character(20), intent(in) :: input

          ns=0
          write(*,*) "open ", input
          open (10, file=input)
          read (10,*) whatever
          read (10,*) command
          Do
             read (10,*,iostat=ios) whatever
             if (ios/=0) exit
             ns=ns+1
          End Do
          close(10)
          write(*,201)  ns
201       format (1x,'electron dynamic steps = ', I10)
       endsubroutine obtain_ns

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

      END PROGRAM TDANALYZE
