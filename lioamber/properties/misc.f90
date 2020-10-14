!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% GET_DEGENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Gets degeneration number (nDeg) for the MO of interest (orb) and the MO      !
! indexes (nDegMO) of those MOs which share degeneration.                      !
subroutine get_degeneration(energ, orb, nDeg, nDegMO)
   ! orb   : The index of the MO we want to get the degeneration of.
   ! energ : The energy of all MOs.
   ! nDeg  : The number of MOs with same energy as "orb". (output)
   ! nDegMO: Indeces of the MOs included in nDeg. (output)
   implicit none
   integer, intent(in)  :: orb
   LIODBLE, intent(in)  :: energ(:)
   integer, intent(out) :: nDeg
   integer, intent(out) :: nDegMO(:)
 
   integer :: iOrb, iCountU, iCountD
   LIODBLE :: criterium, ratio
        
   if (size(energ,1) == 1) then
       nDeg = 1
       nDegMO(1) = 1
       return
   endif

   ratio     = 0.0D0
   nDeg      = 0 
   nDegMO    = 0 
   criterium = 0.00001D0

   ! Assings counting limits in order to avoid exploring all MOs.
   ! Special cases: only one MO, or MO of interest is the first one.
   if (orb > 1) then 
       iCountU = orb
       iCountD = orb - 1
   else
       iCountU = 2
       iCountD = 1
   endif

   ! Performs energy comparisons, upwards and downwards.
   ! The "exit" part is to avoid looping through all functions. If the first
   ! MO we find is not degenerate with the target MO, there is no need to 
   ! continue.
   do iOrb = iCountU, size(energ,1), 1
       ratio = 2.0D0 * &
               abs((energ(iOrb) - energ(orb)) / (energ(iOrb) + energ(orb)))
       if (ratio < criterium) then
           nDeg = nDeg + 1
           nDegMO(nDeg) = iOrb
       else 
           exit
       endif
   enddo
         
   do iOrb = iCountD, 1, -1
       ratio = 2.0D0 * &
               abs((energ(iOrb) - energ(orb)) / (energ(iOrb) + energ(orb)))
       if (ratio < criterium) then 
           nDeg = nDeg + 1
           nDegMO(nDeg) = iOrb
       else
           exit
       endif
   enddo
end subroutine get_degeneration