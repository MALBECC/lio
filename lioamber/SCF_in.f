!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      Subroutine SCF_in(E,qmcoords,clcoords,nsolin,dipxyz)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      use garcha_mod
      REAL*8 , intent(in)  :: qmcoords(3,natom)
      REAL*8 , intent(in)  :: clcoords(4,nsolin)
      nsol=nsolin
      ntatom=nsol+natom

      deallocate (r,v,Em,Rm,pc,nnat)

      allocate (r(ntatom,3),v(ntatom,3),Em(ntatom)
     >,Rm(ntatom),pc(ntatom),nnat(ntatom))
       ngDyn=natom*ng0
       
      if(writexyz) write(18,*) natom
      if(writexyz) write(18,*)

      do i=1,nsol
        n=natom+i
        pc(n)=clcoords(4,i)

        do j=1,3
          r(n,j)=clcoords(j,i)/0.529177D0
        enddo
c        write(18,345) 8,r(n,1),r(n,2),r(n,3)
      enddo


      do i=1,natom
        do j=1,3

c          write(89,*) i, j, qmcoords(i,j)
          r(i,j)= qmcoords(j,i)/0.529177D0
          rqm(i,j)= qmcoords(j,i)/0.529177D0
c          write(87,*) i, j , r(i,j)
        enddo
        write(18,345) Iz(i),qmcoords(:,i)
      enddo

!--------------------------------------------------------------------!
! I am not sure this should be here, but it is the only
! place to put it right now (FFR)
      call lio_init()
!--------------------------------------------------------------------!

      if(OPEN) then 
        call SCFOP(E,dipxyz)
      else
        call SCF(E,dipxyz)
      endif

 345  format(2x,I2,2x,3(f10.6,2x))
      return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
