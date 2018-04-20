!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Steepest Decend algorithm for geometry optimization
!
! Nicolas Foglia, 2017
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

	SUBROUTINE do_steep(E)
	USE garcha_mod, only : Force_cut, Energy_cut, minimzation_steep, n_min_steeps, OPEN, natom, r, rqm, lineal_search, n_points
	IMPLICIT NONE
	real*8, intent(inout) :: E !energy
	real*8 :: Emin, step_size, Fmax, d_E !Energy of previus steep, max displacement (Bohrs) of an atom in each steep, max |force| on an atoms, E(steep i)-E(steep i-1)
	integer :: n !steep number
	integer :: i,j !auxiliars

	double precision, allocatable, dimension(:) :: Energy !array of energy for lineal search
	double precision :: lambda !optimal displacement for minimice energy
	logical :: require_forces ! control force calculations in NO lineal search case
	double precision, dimension(natom, 3) :: r_scrach !positions scratch
	double precision, dimension(3, natom) :: gradient
	logical :: stop_cicle, converged !for stop geometry optimization cicle

	gradient=0.d0

	if ( lineal_search ) write(*,*) "starting Steepest Descend with linear search"
	if ( .not. lineal_search ) write(*,*) "starting Steepest Descend without linear search"

	open (unit=12, file='traj.xyz') ! trajectory file
	open (unit=13, file='optimization.out') !
	write(13,4800)

	n=0
	stop_cicle=.false.
	converged=.false.
	step_size= minimzation_steep
	Fmax=0.d0
	lambda=0.d0
	require_forces=.true.

	call SCF(E)
	Emin=E


	if ( lineal_search ) allocate(Energy(n_points))


	DO WHILE (n .lt. n_min_steeps .and. .not. stop_cicle)
	  write(13,4801) n, E, Fmax, step_size, lambda
	  call save_traj()

	  if ( require_forces) then
	    call dft_get_qm_forces(gradient) !get gradient
	    call max_force(gradient, Fmax) !find max force
	  end if

	  if ( lineal_search ) then
	    call make_E_array(gradient, n_points,Energy, step_size, E,  Fmax) !obtain E(i) with r(i) = r_old - gradient/|gradient| * step_size *i
	    call line_search(n_points,Energy, step_size, lambda ) !predict best lambda that minimice E moving r(lambda) = r_old - gradient/|gradient| * lambda
	  else
	    r_scrach=r
	    lambda=step_size
	  end if

	  call move(lambda, Fmax,gradient)
	  call SCF(E)


	  d_E=abs(E-Emin)
	  if ( lineal_search ) then
	    Emin=E
	    if (lambda .lt. 0.1d0 * step_size) step_size=step_size*0.1d0
	    if (lambda .gt. (dble(n_points)-0.1d0) * step_size) then
	      step_size=step_size*1.5
	      require_forces=.false.
	    else
	      require_forces=.true.
	    end if
	  else
	    if (E .le. Emin) then
	      Emin=E
	      require_forces=.true.
	      step_size=step_size*1.2d0
	    else
	      r=r_scrach
	      require_forces=.false.
	      step_size=step_size*0.5d0
	    end if
	  end if

!convergence criterium
	  if ( step_size .lt. 1D-12) stop_cicle=.true.
	  if ( Fmax .lt. Force_cut .and. d_E .lt. Energy_cut .and. require_forces) then
	    stop_cicle=.true.
	    converged=.true.
	  end if
	  n=n+1
	END DO

        write(13,4801) n, E, Fmax, step_size, lambda
        call save_traj()

	IF (converged) THEN
	  WRITE(*,*) "system converged to required presition in ", n, " steeps"
	ELSE
	  WRITE(*,*) "system has not convergetd to required presition in ", n, " steeps"
	  WRITE(*,*) "Fmax component is: ", Fmax, "and energy diferences is: ", d_E
	END IF

	close(12)
	close(13)

	if ( lineal_search ) deallocate(Energy)
 4800   FORMAT(6x, "steep", 7x, "energy", 14x, "Fmax ", 9x, "steep size", 10x, "lambda")
 4801   FORMAT(2x,"!",2x, i4, 2x, 4(f16.10,2x))
	END SUBROUTINE do_steep


	subroutine max_force(gradient, Fmax)
	use garcha_mod, only : natom
	implicit none
	double precision, intent(out) :: Fmax
	double precision, intent(in) :: gradient(3,natom)
	double precision :: F_i
	integer :: i
	Fmax=0.d0
        do i=1, natom
          F_i=gradient(1,i)**2 + gradient(2,i)**2 + gradient(3,i)**2
          F_i=sqrt(F_i)
          if (F_i .gt. Fmax) Fmax=F_i
        end do
	return
	end subroutine max_force


	subroutine make_E_array(gradient, n_points,Energy, step_size, E, Fmax)
!generate E(i) moving system r_new=r_old - i*step_size*gradient
	use garcha_mod, only : r, natom, rqm
   use fileio_data, only: verbose
	implicit none
	double precision, intent(in) :: gradient(3,natom), step_size, Fmax
        double precision, intent(inout) :: E
	integer, intent(in) :: n_points
	double precision, intent(out) :: Energy(n_points)
	double precision, dimension(natom, 3) :: r_ini
	double precision :: a
	double precision :: max_move
	integer :: i,j,k

!define step that normalice max force
	max_move=step_size/Fmax

!generates Energy(i)
	r_ini=r
	Energy(1)=E
	do i=2, n_points
	  do j=1, natom  !move positions along gradient direction
	    do k=1, 3
	      r(j,k)=r_ini(j,k)-max_move*dble(i)*gradient(k,j)
	    end do
	  end do
	  rqm=r
	  call SCF(E)
	  Energy(i)=E
	end do

	if(verbose .gt. 1) then !ponerle luego un nivel alto de verbose
	  do i=1, n_points
	    write(13,5008) Energy(i)
	  end do
	end if

	r=r_ini
	rqm=r
	return
 5008 FORMAT(2x,"Lineal search Energy ",2x,f16.10)
	end subroutine make_E_array



        subroutine line_search(n_points,Energy, step_size, lambda )
!if minimun value of E is E(n_points) return lambda = step_size * n_points
!if minimun value of E is E(1) return lambda = 0.d0
!else returns an expected lambda that minimice E(lambda) using a parabolic interpolation
        implicit none
        integer, intent(in) :: n_points
        double precision, dimension(n_points), intent(in) :: Energy
        double precision, intent(in) :: step_size
        double precision, intent(out) :: lambda
        integer :: i
        integer :: min_Energy_position
        double precision :: dE1, dE2, modif_fac

        if (n_points .le. 2) then
	  write(*,*) "wrong n_points in lineal search, n_points need to be > 2"
	  stop
	end if
!find min value in Energy elements
        min_Energy_position=1
        do i=2, n_points
           if (Energy(i) .lt. Energy(min_Energy_position)) min_Energy_position=i
        end do

        if (min_Energy_position .eq. 1) then
          lambda=0.d0
          return
        elseif (min_Energy_position .eq. n_points) then
          lambda=step_size*dble(n_points)
          return
        end if

        dE2=abs(Energy(min_Energy_position) - Energy(min_Energy_position+1))
        dE1=abs(Energy(min_Energy_position) - Energy(min_Energy_position-1))

        modif_fac=step_size*(dE2-dE1)/(dE1+dE2)
        lambda=step_size*dble(min_Energy_position) - 0.5d0 * modif_fac
        return
        end subroutine line_search


        SUBROUTINE move(lambda, Fmax, gradient) !calculate new positions
	USE garcha_mod, only : r, rqm,  natom
        IMPLICIT NONE
        INTEGER :: i,j
        DOUBLE PRECISION, INTENT(IN) :: lambda, Fmax, gradient(3,natom)
	DOUBLE PRECISION :: a
        a=lambda/Fmax
        DO i=1,natom
          DO j=1,3
             r(i,j)= r(i,j)-a*gradient(j,i)
          END DO
        END DO
	rqm=r
        END SUBROUTINE move

	SUBROUTINE save_traj()
	USE garcha_mod, only : IZ, r, natom
	IMPLICIT NONE
	integer :: i
	REAL*8 :: b_2_ang
	b_2_ang=0.52917725D0
        write(12,*) NATOM
        write(12,*)
        DO i=1,NATOM
          write(12,5000) IZ(i), r(i,1)*b_2_ang, r(i,2)*b_2_ang, r(i,3)*b_2_ang
        END DO
 5000 FORMAT(2x,i2,2x,3(f16.10,2x))
	END SUBROUTINE save_traj
