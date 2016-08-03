      SUBROUTINE SEEK_NaN(VEC_TO_TEST,INICIO, FIN,FRASE)
!subroutine for seek NaN in vectors
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: VEC_TO_TEST(*) !vector to analize
        INTEGER, INTENT(IN) :: INICIO, FIN  !range of analize in vector
        CHARACTER (LEN=*), INTENT(IN) :: FRASE !output phrase in case of NaN
        INTEGER :: iNick

        IF (INICIO .GT. FIN) THEN
          WRITE(*,*) "Error: Inicio mayor a Fin en test"
          WRITE(*,*) FRASE
          STOP
        END IF

        do iNick=INICIO,FIN
         if (VEC_TO_TEST(iNick) .ne. VEC_TO_TEST(iNick)) then
          write(*,*) "NAN en ",FRASE, iNick
          stop
         end if
        end do
      ENDSUBROUTINE SEEK_NaN

