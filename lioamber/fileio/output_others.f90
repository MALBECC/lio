!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% WRITE_OUTPUT.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains several output-file printing routines. Currently includes:!
! * write_dipole     (handles dipole moment printing)                          !
! * write_forces     (handles grandient printing to output)                    !
! * write_fukui      (handles Fukui function printing to output)               !
! * write_orbitals   (prints orbitals and energies to output)                  !
! * write_population (handles population/charge printing to output)            !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% WRITE_DIPOLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Sets variables up and calls dipole calculation.                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_dipole(dipxyz, u, uid, header)
    use fileio_data, only : style

    implicit none
    real*8 , intent(in) :: dipxyz(3), u
    integer, intent(in) :: uid
    logical, intent(in) :: header


    open(unit = uid, file = "dipole_moment")
    if (style) then
        if (header) then
            write(UID,8698)
            write(UID,8699)
            write(UID,8700)
            write(UID,8701)
            write(UID,8702)
        else
            write(UID,8704) dipxyz(1), dipxyz(2), dipxyz(3), u
        endif
    else
        if (header) then
            write(UID,*)
            write(UID,*) '#DIPOLE MOMENT, X Y Z COMPONENTS AND NORM (DEBYES)'
            write(UID,*)
        else
            write(UID,*) dipxyz(1), dipxyz(2), dipxyz(3), u
        endif
    endif

 8698 FORMAT(4x,"╔════════════════",&
      "═════════════════════",&
      "═════════════════════","═════╗")

 8699 FORMAT(4x,"║                         Dipole Moment            ", &
      "             ║")
 8700 FORMAT(4x,"╠═══════════════╦", &
      "═══════════════╦═════",       &
      "══════════╦══════════",       &
      "═════╣")
 8701 FORMAT(4x,"║       ux      ║       uy      ║       uz     ",&
      " ║       u       ║")
 8702 FORMAT(4x,"╠═══════════════╬", &
      "═══════════════╬═════",       &
      "══════════╬══════════",       &
      "═════╣")
 8704 FORMAT(4x,4("║"F13.9,2x),"║")

    return
end subroutine write_dipole
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% WRITE_DIPOLE_TD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Sets variables up and calls dipole calculation.                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_dipole_td(dipxyz, time, uid)
    implicit none
    real*8 , intent(in) :: dipxyz(3), time
    integer, intent(in) :: uid

    write(UID,100) time, dipxyz(1), dipxyz(2), dipxyz(3)
100 format (e15.8,' ', e15.8,' ',e15.8,' ',e15.8)
    return
end subroutine write_dipole_td
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% WRITE_DIPOLE_TD_HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Sets variables up and calls dipole calculation.                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_dipole_td_header(time_step, fx, fy, fz, uid)
    implicit none
    real*8 , intent(in) :: time_step, fx, fy, fz
    integer, intent(in) :: uid

    write(UID, 100) time_step, fx, fy, fz
100 format ('# ',e12.5,' ', e12.5,' ',e12.5,' ',e12.5)
    return
end subroutine write_dipole_td_header

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% WRITE_FORCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Calculates forces for QM and MM regions and writes them to output.           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_forces(dxyz, natom, offset, uid)
    implicit none
    integer, intent(in) :: uid, natom, offset
    real*8 , intent(in) :: dxyz(3, natom+offset)
    integer             :: kk
    do kk=offset+1, offset+natom
        write(uid,100) kk, dxyz(1, kk), dxyz(2, kk), dxyz(3, kk)
    enddo

100 format (I5,2x,f10.6,2x,f10.6,2x,f10.6)

    return
end subroutine write_forces
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% WRITE_FUKUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Writes Fukui function and local softness to output.                          !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_fukui(fukuiNeg, fukuiPos, fukuiRad, N, Iz, soft)

    implicit none
    integer, intent(in) :: N, Iz(N)
    real*8 , intent(in) :: fukuiNeg(N), fukuiPos(N), fukuiRad(N), soft
    integer :: i

    write(*,*) "Global Softness (A.U.):  ", soft
    write(*,*) "N", "Fukui-", "Fukui+", "Fukui0", "Local Softness (A.U.)"
    do i=1, N
        write(*,*) Iz(i), fukuiNeg(i), fukuiPos(i), fukuiRad(i), &
                   abs(soft*fukuiRad(i))
    enddo

end subroutine write_fukui
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% WRITE_ORBITALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Prints orbital energies and coefficients.                                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_orbitals(UID)
    use garcha_mod, only : NCO, M, X, Eorbs

    implicit none
    integer, intent(in) :: UID
    integer             :: nn, ll

    write(UID,*) 'ORBITAL COEFFICIENTS AND ENERGIES, CLOSED SHELL'
    do nn = 1, NCO
       write(UID, 850) nn, Eorbs(nn)
       write(UID, 400) (X(ll, 2*M +nn), ll = 1, M)
    enddo
    do nn = NCO+1, M
       write(UID, 851) nn, Eorbs(nn)
       write(UID, 400) (X(ll, 2*M +nn), ll = 1, M)
    enddo

400 format(4(E14.7E2, 2x))
850 format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7)
851 format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7, '(NON OCC.)')

    return
endsubroutine write_orbitals
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% WRITE_POPULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Writes Mulliken/Löwdin charges to output.                                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_population(UID, N, q0, q, pop)

    implicit none
    integer, intent(in) :: UID, N, q0(N), pop
    real*8 , intent(in) :: q(N)

    real*8  :: qtotal
    integer :: i

    call g2g_timer_start('mulliken_write')
    qtotal=0.d0

    write(UID,*)
    write(UID,300)
    if (pop.eq.0) write(UID,301)
    if (pop.eq.1) write(UID,309)
    write(UID,302)
    write(UID,303)
    write(UID,304)
    do i=1,N
        qtotal=qtotal+q(i)
        write(UID,305) i, q0(i), q(i)
    enddo
    write(UID,306)
    write(UID,307) qtotal
    write(UID,308)
    write(UID,*)

    call g2g_timer_stop('mulliken_write')

300 FORMAT(8x,"╔════════════════", &
    "═════════════════╗")
301 FORMAT(8x,"║   MULLIKEN POPULATION ANALYSIS  ║")
302 FORMAT(8x,"╠════════╦═══════════╦════════════╣")
303 FORMAT(8x,"║ ATOM # ║ ATOM TYPE ║ POPULATION ║")
304 FORMAT(8x,"╠════════╬═══════════╬════════════╣")
305 FORMAT(8x,"║",2x,i3,3x,"║"3x,i3,5x,"║",1x,F10.7,1x,"║")
306 FORMAT(8x,"╚════════╬═══════════╬════════════╣")
307 FORMAT(8x,"         ║   TOTAL   ║",1x,F10.7,1x,"║")
308 FORMAT(8x,"         ╚═══════════╩════════════╝")
309 FORMAT(8x,"║    LÖWDIN POPULATION ANALYSIS   ║")
    return
end subroutine write_population
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
