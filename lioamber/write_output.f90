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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% WRITE_RESTART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%m
