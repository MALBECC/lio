subroutine energy_specific(eigvec,eigval,Sdim,Nstat)
! This routine controls the Energy Specific TDA
! Paper: DOI= 10.1021/ct200485x
! energy_min = lowest energy in LR solution

use excited_data, only: estda, energy_min
   implicit none

   integer, intent(in) :: Sdim, Nstat
   LIODBLE, intent(inout) :: eigvec(Sdim,Sdim), eigval(Sdim)

   integer :: ii, start, cont
   LIODBLE, dimension(:,:), allocatable :: vec_temp
   LIODBLE, dimension(:)  , allocatable :: val_temp

   if ( .not. estda ) return

   allocate(vec_temp(Sdim,Sdim),val_temp(Sdim))
   vec_temp = eigvec; val_temp = eigval

   start = 0
   do ii=1,Sdim
      if (eigval(ii) > energy_min) then
         start = ii
         exit
      endif
   enddo
   if ( start == 0 ) then
      print*, "There are not eigenvalues lower than threshold"
      stop
   endif

   if ( (Sdim-start) < Nstat ) then
      print*, "There are not enought eigenvectors with the energy specific"
      stop
   endif

   eigval = 0.0d0; eigvec = 0.0d0
   cont = 1
   do ii=start, start+Nstat
      eigval(cont)   = val_temp(ii)
      eigvec(:,cont) = vec_temp(:,ii)
      cont = cont + 1
   enddo

   deallocate(val_temp,vec_temp)
end subroutine energy_specific

subroutine estda_window(deltaE,ind,Ndim,start,vecnum,Sdim,maxIt)
use excited_data, only: energy_min, d_energy, window
   implicit none
   
   integer, intent(in) :: Ndim
   integer, intent(in) :: ind(Ndim)
   LIODBLE, intent(in) :: deltaE(Ndim)
   integer, intent(inout) :: start, vecnum, Sdim, maxIt

   logical :: find
   integer :: ii, nwin, idwin
   LIODBLE :: oldE, difE, Ewin, MINDIF
   LIODBLE, dimension(:), allocatable :: Estartw

   print*, " "
   print*, "Energy Specific TDA"
   write(*, "(1X,A,F8.4)") "energy_min= ",energy_min

   ! How many windos has the system
   oldE = deltaE(ind(1)); nwin = 0
   do ii=2,Ndim
      difE = deltaE(ind(ii)) - oldE
      if ( difE > d_energy ) nwin = nwin + 1
      oldE = deltaE(ind(ii))
   enddo
   write(*, "(1X,A,F8.4,A,I3,A)") "With threshold ", d_energy, " You have ", nwin, " windows"

   ! Get start Energy window
   allocate(Estartw(nwin)); Estartw = 0.0d0
   oldE = deltaE(ind(1)); nwin = 0
   do ii=2,Ndim
      difE = deltaE(ind(ii)) - oldE
      if ( difE > d_energy ) then
         nwin = nwin + 1
         Estartw(nwin) = deltaE(ind(ii))
      endif
      oldE = deltaE(ind(ii))
   enddo
   write(*,"(1x,A)") "Start Energy Window [Ha.]"
   do ii=1,nwin 
      write(*,"(1X,F12.6)",ADVANCE="NO") Estartw(ii)
   enddo
   write(*,*)

   ! Select what window will be calculating
   MINDIF = dabs(energy_min-Estartw(1)) + 1.0d0
   idwin  = 0
   Ewin   = 0.0d0
   do ii=1,nwin
      difE = dabs(energy_min-Estartw(ii))
      if ( difE < MINDIF ) then
         MINDIF = difE
         idwin  = ii
         Ewin   = Estartw(ii)
      endif
   enddo
   if ( window /= -1 ) then
      idwin = window
      Ewin  = Estartw(idwin)
   endif
   write(*,"(1X,A,I3,A,F12.8,A)") "We'll calculate the window ", idwin, ", ", Ewin, " Ha."
   energy_min = Ewin
   d_energy   = -d_energy/20.d0

   ! Generating Trials vectors in this window
    do ii=1,Ndim
       if ( deltaE(ind(ii))>(energy_min+d_energy) ) then
          start = ii
          find = .true.
          exit
       endif
    enddo
    if ( .not. find ) then
       print*, "There is not energy root above the threshold"
       stop
    endif

    ! Check new dimensions in ES-TDA
    if ( (start+vecnum) >= Ndim ) then
       vecnum=Ndim-start
       Sdim  =vecnum
       maxIt = 1
    endif
    print*, " "

    deallocate(Estartw)
end subroutine estda_window
