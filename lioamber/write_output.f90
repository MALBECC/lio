!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% PRINT_ORBITALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Prints orbital energies and coefficients.                                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine print_orbitals(UID)
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
endsubroutine print_orbitals
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% WRITE_DIPOLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Sets variables up and calls dipole calculation.                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_dipole(dipxyz, u, uid)
    use garcha_mod, only : style

    implicit none
    integer, intent(in) :: uid
    real*8 , intent(in) :: dipxyz(3), u

    open(unit = uid, file = "dipole_moment")
    if (style) then
        write(UID,8698)
        write(UID,8699)
        write(UID,8700)
        write(UID,8701)
        write(UID,8702)
        write(UID,8704) dipxyz(1), dipxyz(2), dipxyz(3), u
    else
        write(UID,*)
        write(UID,*) 'DIPOLE MOMENT, X Y Z COMPONENTS AND NORM (DEBYES)'
        write(UID,*)
        write(UID,*) dipxyz(1), dipxyz(2), dipxyz(3), u
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
!%% WRITE_FORCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Calculates forces for QM and MM regions and writes them to output.           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_forces(dxyz, natom, offset, uid)
    implicit none
    integer, intent(in) :: uid, natom, offset
    real*8 , intent(in) :: dxyz(natom, 3)
    integer             :: kk

    do kk=offset, offset+natom
        write(uid,100) kk, dxyz(kk, 1), dxyz(kk, 2), dxyz(kk, 3)
    enddo

100 format (I5,2x,f10.6,2x,f10.6,2x,f10.6)

    return
end subroutine write_forces
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% POPULATION_WRITE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Writes Mulliken/Löwdin charges to output.                                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine population_write(UID, N, q0, q, pop)

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

300 FORMAT(8x,"╔═════════════════&
    ════════════════╗")
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
end subroutine population_write
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% WRITE_RESTART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Writes a coefficient restart file.                                           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_restart(UID)
    use garcha_mod, only : M, X, NCO, indexii
 
    implicit none
    integer, intent(in) :: UID
    integer             :: ll, nn

    call g2g_timer_sum_start('restart write')
    rewind UID
    do ll = 1, M
        do nn = 1, M
            X(indexii(ll), M +nn) = X(ll, 2*M +nn)
        enddo
    enddo

    do ll = 1, M
         write(UID,400) (X(ll, M +nn), nn = 1, NCO)
    enddo

    call g2g_timer_sum_stop('restart write')

400 format(4(E14.7E2, 2x))
    return
end subroutine write_restart
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
