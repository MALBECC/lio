!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine get_restrain_forces(dxyzqm, f_r)
!--------------------------------------------------------------------!
!This routine calculate forces for a distance restrain
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       use garcha_mod,only:natom,r, number_index, number_restr,         &
       restr_index, restr_pairs, restr_k, restr_r0, restr_w
       use fileio_data, only: verbose
       implicit none
       integer :: i,j, k, l !auxiliar
       real*8, dimension(natom,natom) :: distx,disty,distz
       real*8,intent(inout) :: dxyzqm(3,natom)
       real*8, intent(out) :: f_r
       real*8 :: k_force, r_eq, W_distance, distance
       real*8 :: Fx, Fy, Fz
       integer :: ai, aj
!--------------------------------------------------------------------!

!Calculate distances
        DO i=1,natom
          DO j=1,natom
            distx(i,j)=r(i,1)-r(j,1)
            disty(i,j)=r(i,2)-r(j,2)
            distz(i,j)=r(i,3)-r(j,3)
          END DO
        END DO

	IF ( verbose.gt.0) THEN
         open(unit=1001,file="lio.restrain.out")
         write(1001,*) "  index        Coord. Value               Force"
        END IF

        DO l=0, number_index-1 ! Loop over indexes, for distance combination
          W_distance=0.d0
	  k_force=0.d0
!first calculate prefactor
          DO k=1, number_restr !Loop over restrains
            IF (restr_index(k) .EQ. l) THEN !restranis that will treath together
              ai=restr_pairs(1,k) !number of first atom fo restrain
              aj=restr_pairs(2,k) !number of second atom of restrain
              distance= distx(ai,aj)**2 + disty(ai,aj)**2 + distz(ai,aj)**2
              distance= distance**0.5
              W_distance=W_distance + distance*restr_w(k)
              k_force=restr_k(k)
              r_eq=restr_r0(k)
            END IF
          END DO

	  f_r=-k_force*(W_distance-r_eq)

        IF ( verbose.gt.0) THEN
         write(1001,5500) l, W_distance, f_r
        END IF



          DO k=1, number_restr !Loop over restrains
            IF (restr_index(k) .EQ. l) THEN
              ai=restr_pairs(1,k) !number of first atom of restrain
              aj=restr_pairs(2,k) !number of second atom of restrain
              distance= distx(ai,aj)**2 + disty(ai,aj)**2 + distz(ai,aj)**2
              distance= distance**0.5

	      IF (distance .eq. 0) STOP "Distance is 0 for 2 atoms"

	      Fx=f_r*restr_w(k)*distx(ai,aj)/distance
	      Fy=f_r*restr_w(k)*disty(ai,aj)/distance
	      Fz=f_r*restr_w(k)*distz(ai,aj)/distance


!signs inverted here because dxyzqm stores gradients
	      dxyzqm(1,ai)=dxyzqm(1,ai)-Fx
              dxyzqm(2,ai)=dxyzqm(2,ai)-Fy
              dxyzqm(3,ai)=dxyzqm(3,ai)-Fz

              dxyzqm(1,aj)=dxyzqm(1,aj)+Fx
              dxyzqm(2,aj)=dxyzqm(2,aj)+Fy
              dxyzqm(3,aj)=dxyzqm(3,aj)+Fz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            END IF
          END DO
        END DO

        IF ( verbose.gt.0) THEN
         close (1001)
        END IF

 5500  format(4x,I2,2x,2(f22.16,2x))
       end subroutine get_restrain_forces


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine get_restrain_energy(E_restrain)
!--------------------------------------------------------------------!
!This routine calculate energy for a distance restrain
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       use garcha_mod,only:natom,r, number_index, number_restr,      &
       restr_index, restr_pairs, restr_k, restr_r0, restr_w
       implicit none
       integer :: i,j,k,l !auxiliar
       real*8, dimension(natom,natom) :: distx,disty,distz
       real*8,intent(out) :: E_restrain
       real*8 :: k_force, r_eq, W_distance, distance
       integer :: ai, aj
!--------------------------------------------------------------------!
       E_restrain=0.d0

!Calculate distances
        DO i=1,natom
          DO j=1,natom
            distx(i,j)=r(i,1)-r(j,1)
            disty(i,j)=r(i,2)-r(j,2)
            distz(i,j)=r(i,3)-r(j,3)
          END DO
        END DO

	DO l=0, number_index ! Loop over indexes, for distance combination
	  W_distance=0.d0
	  DO k=1, number_restr !Loop over restrains
	    IF (restr_index(k) .EQ. l) THEN
	      ai=restr_pairs(1,k) !number of first atom fo restrain
	      aj=restr_pairs(2,k) !number of second atom of restrain
	      distance= distx(ai,aj)**2 + disty(ai,aj)**2 + distz(ai,aj)**2
	      distance= distance**0.5
	      W_distance=W_distance + distance*restr_w(k)
	      k_force=restr_k(k)
	      r_eq=restr_r0(k)
	    END IF
	  END DO
	  E_restrain=E_restrain + k_force/2.d0 *(W_distance-r_eq)**2
	  k_force=0.d0
	END DO
	RETURN
       end subroutine get_restrain_energy



	SUBROUTINE read_restrain_params()
	USE garcha_mod, ONLY : number_restr, restr_pairs, restr_index, restr_k, restr_w, restr_r0, number_index
	IMPLICIT NONE
        LOGICAL :: exist_restr_file, add
	INTEGER :: i, j
        INQUIRE(FILE="lio.restrain", EXIST=exist_restr_file)
        IF ( .NOT. exist_restr_file) STOP "lio.restrain file not found" !check existence of restrain file

        OPEN(UNIT=25,FILE="lio.restrain")
	READ(25,*)
	DO i=1, number_restr
	READ(25,*) restr_pairs(1,i), restr_pairs(2,i), restr_index(i), restr_k(i), restr_w(i), restr_r0(i)
	END DO
	CLOSE(25)

	!calcula cuantos index hay
	number_index=0
	DO j=0, 10
	  add=.FALSE.
	  DO i=1, number_restr
	    IF (j.EQ.restr_index(i)) add=.TRUE.
	  END DO
	  IF (add) number_index = number_index +1
	END DO

	RETURN
	END SUBROUTINE read_restrain_params
