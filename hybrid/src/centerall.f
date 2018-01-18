c----------------------------------------------------------------------------
c    Centers quantum subsystem inside the cell
c----------------------------------------------------------------------------
c
      subroutine centermol(na_u,xa,rclas,ucell,natot)
                       
      implicit none

      integer i,j,na_u,natot
      double precision rclas(3,natot), xa(3,na_u), ucell(3,3),
     . x(3), d0(3), min(3), max(3), delta(3)

      d0=0.0D0
 
      min(1:3)=xa(1:3,1)
      max(1:3)=xa(1:3,1)
     
      do j=1,3 
       do i=1,na_u
        if (xa(j,i).gt.max(j)) max(j)=xa(j,i)
        if (xa(j,i).lt.min(j)) min(j)=xa(j,i)
       enddo
      enddo

      do i=1,3
       delta(i)=max(i)-min(i)
       d0(i)=(ucell(i,i)-delta(i))/2.
      enddo

      do i=1,na_u
       xa(1:3,i)=xa(1:3,i)-min(1:3)+d0(1:3)
      enddo
        
      do i=na_u+1,natot
       rclas(1:3,i)=rclas(1:3,i)-min(1:3)+d0(1:3)
      enddo

      do i=1,na_u
       rclas(1:3,i)=xa(1:3,i)
      enddo

      write(6,'(/a)') 'center: System centered in cell'

      return
      end

c----------------------------------------------------------------------------
c    Centers quantum subsystem inside the cell 
c    when dynamics drives it near a border
c----------------------------------------------------------------------------
c
      subroutine centerdyn(na_u,rclas,ucell,natot)
                       
      implicit none

      integer i,j,na_u,natot
      double precision rclas(3,natot), ucell(3,3),
     . x(3), d0(3), min(3), max(3), delta(3), lbord, dbord(3)

      logical ctr
      lbord = 2.0

      d0=0.0D0
      dbord=0.0D0
      ctr=.false.
 
      min(1:3)=rclas(1:3,1)
      max(1:3)=rclas(1:3,1)
     
      do j=1,3 
       do i=1,na_u
        if (rclas(j,i).gt.max(j)) max(j)=rclas(j,i)
        if (rclas(j,i).lt.min(j)) min(j)=rclas(j,i)
       enddo
      enddo

       do i=1,3
        dbord(i)=ucell(i,i)-max(i)
        if (min(i).lt.lbord.or.dbord(i).lt.lbord) ctr=.true.
       enddo

       if (ctr) then

       write(6,'(/a)') 'center: System centered in cell'

       do i=1,3
        delta(i)=max(i)-min(i)
        d0(i)=(ucell(i,i)-delta(i))/2.
        if (d0(i).lt.lbord) then
         write(6,'(/a)') 'center: QM atoms too close to cell border'
        endif
       enddo

       do i=1,natot
        rclas(1:3,i)=rclas(1:3,i)-min(1:3)+d0(1:3)
       enddo
 
       endif 

       return
       end

c----------------------------------------------------------------------------
