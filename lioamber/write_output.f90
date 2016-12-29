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
