
      subroutine jar(f) 
          use garcha_mod
c
      implicit real*8 (a-h,o-z)
      dimension f(natom,3)
c      dimension jatom(2,100),coef(100)i,dist(100,3),distt(100)
      open(unit=33,file='jarz')

      read(33,*)
      read(33,*) ndis
       if (ndis.gt.99) then
       write(*,*) 'ndist > 100'
       stop
      endif
      read(33,*)
      do i=1,ndis
      read(33,*) jatom(1,ndis), jatom(2,ndis)
      enddo
      read(33,*)
      do i=1,ndis
      read(33,*) coef(ndis)
      enddo
      read(33,*)
      read(33,*)  Kjar
      read(33,*)
      read(33,*)  xini,xfinal,steps
      polin=0.
      xjar=(xfinal-xini)/steps*(npas-1)+xini
      do i=1,ndist
      do j=1,3
        dist(i,j)= (r(j,jatom(1,i))-r(j,jatom(2,i)))

         enddo
        distt(i)=sqrt(dist(i,1)**2+dist(i,2)**2+dist(i,3)**2)

       polin=coef(i)*distt(i)
      enddo

      factor=kjar*(polin-xjar)

      do i=1,ndis
       do j=1,3

        f(jatom(1,ndis),j)=f(jatom(1,ndis),j)-factor*dist(i,j)/distt(i)
        f(jatom(2,ndis),j)=f(jatom(2,ndis),j)+factor*dist(i,j)/distt(i)
        
       enddo
      enddo


       return
       end
