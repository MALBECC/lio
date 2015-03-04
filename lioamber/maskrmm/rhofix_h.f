!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!  RhoFix // RhoMess
!------------------------------------------------------------------------------!
!
! In the RMM vector, Rho is stored in packed storage but
! each non diagonal position has dobule the real value.
! So before putting a density matrix in RMM, non diagonal
! positions need to be x2 (messrho) and when taking it
! out of RMM, the non diagonal positions need to be
! divided by two (fixrho).
!
! Note that fixrho modifies the matrix, whereas mess rho
! modifies the vector.
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       interface rhofix
         module procedure rhofix_mat
         module procedure rhofix_vec
       end interface
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
