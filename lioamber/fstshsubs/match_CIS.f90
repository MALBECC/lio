subroutine match_CIS(Wnow,Wold,Esup,ndets,all_states,current_state)
use fstsh_data, only: tsh_file, tsh_nucStep, coef_Stat
   implicit none

   integer, intent(in)    :: ndets, all_states
   integer, intent(inout) :: current_state
   LIODBLE, intent(in)    :: Esup(all_states)
   LIODBLE, intent(in)    :: Wnow(ndets,all_states), Wold(ndets,all_states)

   integer :: ii, jj
   LIODBLE :: norm, deltaE
   TDCOMPLEX :: cc
   LIODBLE, allocatable :: vec_norm(:)

   if ( tsh_nucStep == 0 ) return

   allocate(vec_norm(all_states))

   do ii=1,all_states
      norm = 0.0d0
      do jj=1,ndets
         norm = norm + abs( Wnow(jj,current_state) * Wold(jj,ii) )
      enddo
      vec_norm(ii) = norm
      write(tsh_file,"(2X,F10.5)", ADVANCE="NO") norm
   enddo
   write(tsh_file,*) " "

   do ii=1,all_states
      deltaE = abs( Esup(current_state) - Esup(ii) )           ! 0.0048 ha. == 3 Kcal / mol
      if ( vec_norm(current_state) < vec_norm(ii) .and. deltaE < 0.0048d0 ) then
         write(tsh_file,*) "Switch Coef", current_state, "->", ii
         cc = coef_Stat(1,current_state); coef_Stat(1,current_state) = coef_Stat(1,ii)
         coef_Stat(1,ii) = cc
         !current_state = ii
      endif
   enddo
      
   deallocate(vec_norm)
end subroutine match_CIS
