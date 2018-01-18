c Reads Link Atoms parameters
	subroutine link1(numlink,linkat,linkqm,linkmm,
     .  linkqmtype,linkmm2,ng1,namber,na,qmtype,r)
c namber atomos clasicos
c na atomos cuanticos

	use fdf
	use sys
	implicit none
	integer na,namber,i,j,k,m,iunit,numlink,linkqm(15,4),linkmm(15,4),
     .  linkat(15),ng1(namber,6),linkmm2(15,4,3),min(na)
	character exp,linkqmtype(15,4)*4,qmtype(na)*4,ch4*4,ch1*1,ch2*2
	double precision r(3,na+namber),rmin(na),r1,dist


	linkqm=0
	linkmm=0
	linkmm2=0

c check number of linkatoms
	if(numlink.gt.15) then
	  call die('linkatom: number of link atoms must not exeed 15')
!no se por que pusieron este limite, pero hayq  hacerlo mas flexible en el futuro, Nick
	endif

c assignation of linkatoms
	do i=1,numlink
	  linkat(i)=na-numlink+i
	enddo

c assignation of CQM and CMM 
	do i=1,na
	  rmin(i)=dist(r(1,i),r(2,i),r(3,i),r(1,na+1),r(2,na+1),r(3,na+1))
	  min(i)=na+1
	  do j=na+2,na+namber
	    r1=dist(r(1,i),r(2,i),r(3,i),r(1,j),r(2,j),r(3,j))
	    if(r1.le.rmin(i)) then
	      rmin(i)=r1
	      min(i)=j
	    endif
	  enddo
	enddo

	do i=1,na
	  ch4=qmtype(i)
	  ch1=ch4(1:1)
	  if(ch1.ne.'C'.and.ch1.ne.'N') rmin(i)=20
!solo pone link atoms a C y N cuanticos
	enddo


	do k=1,numlink
	  linkqm(k,1)=1
	  do i=2,na
	    if(rmin(i).le.rmin(linkqm(k,1))) then 
	      linkqm(k,1)=i
	    endif
	  enddo
	  if(rmin(linkqm(k,1)).gt.4) then
	    write(6,*) k,linkqm(k,1),rmin(linkqm(k,1))
	    write(6,*) 'Wrong Link Atom CQM atom....Check geometry'
	    STOP
!el atomo QM esta muy lejos de un atomo MM
	  endif
	  rmin(linkqm(k,1))=20
	  linkmm(k,1)=min(linkqm(k,1))-na
	enddo

c assignation of QM 1st neighbors
	rmin=0.0
	min=0
	do i=1,numlink
	m=0
c sp2 carbon
	ch4=qmtype(linkqm(i,1))
        ch1=ch4(1:1)
        if(ch1.eq.'C'.or.ch1.eq.'N') m=3
c sp3 carbon
        ch4=qmtype(linkqm(i,1))
        ch2=ch4(1:2)
	if(ch2.eq.'CT') m=4
	if(m.eq.0) then
	write(6,*) 'Wrong LA QM atom type....Check parameters'
	STOP
	endif
c loop over CQM neighbors
	  do k=2,m
		rmin(i)=10
		min(i)=0
		do j=1,na-numlink
	if(j.eq.linkqm(i,1).or.j.eq.linkqm(i,2).or.j.eq.linkqm(i,3)) then
	goto 5
	endif
		r1=dist(r(1,linkqm(i,1)),r(2,linkqm(i,1)),r(3,linkqm(i,1)),
     .          r(1,j),r(2,j),r(3,j))      
		if(r1.le.rmin(i)) then
		rmin(i)=r1
		min(i)=j
		endif
 5		enddo
	   if(rmin(i).gt.4) then
	   write(6,*) 'Wrong LA QM neighbor....Check geometry'
	   STOP 
	   endif
	    linkqm(i,k)=min(i)
	  enddo
	enddo

c checking CQM neighbors
        do i=1,numlink
c sp2 carbon
        ch4=qmtype(linkqm(i,1))
        ch1=ch4(1:1)
        if(ch1.eq.'C') then
      if(linkqm(i,2).eq.0.and.linkqm(i,3).eq.0.or.linkqm(i,2).eq.0.and.
     .linkqm(i,4).eq.0.or.linkqm(i,3).eq.0.and.linkqm(i,4).eq.0) then
	write(6,*) 'Wrong QM neighbor number for a sp2 carbon'   
	STOP
	endif
	endif
c sp3 carbon
        ch4=qmtype(linkqm(i,1))
        ch2=ch4(1:2)
        if(ch2.eq.'CT') then
	if(linkqm(i,2).eq.0.or.linkqm(i,3).eq.0.or.linkqm(i,4).eq.0) then    
        write(6,*) 'Wrong QM neighbor number for a sp3 carbon' 
	STOP
	endif
	endif
	enddo
	
c asignation of QM link atoms types
        do i=1,numlink
        do j=1,4
	if(linkqm(i,j).eq.0) then
	   linkqmtype(i,j)='XX'
	else
	  linkqmtype(i,j)=qmtype(linkqm(i,j))
	end if

	enddo
	enddo

c assignation of MM 1st neighbors
        do i=1,numlink
          do k=2,4
          linkmm(i,k)=ng1(linkmm(i,1),k-1)
          enddo
        enddo	

c assignation of MM 2nd neighbors
        do i=1,numlink
        do j=2,4
	if(linkmm(i,j).eq.0) goto 10
        m=1
        do k=1,4
        if(ng1(linkmm(i,j),k).ne.linkmm(i,1)) then
        linkmm2(i,j,m)=ng1(linkmm(i,j),k)
        m=m+1
        endif
        enddo
 10     enddo
        enddo	

c write Link Atoms params in file
        write(6,'(/,a)') 'hybrid: Link atom parameters:'
	do i=1,numlink
        write(6,'(a,2x,1I6)') 'linkatom:',linkat(i)
        write(6,'(a,2x,4I6)') 'qmatoms: ',(linkqm(i,j),j=1,4)
C       write(6,'(a,2x,4A4)') 'qmtypes: ',(linkqmtype(i,j),j=1,4)
        write(6,'(a,2x,4I6)') 'mmatoms: ',(linkmm(i,j),j=1,4)
C       write(6,'(a,2x,9I6)') 'mmneihg: ',((linkmm2(i,j,k),k=1,3),j=2,4)
	enddo

	end

c*************************************************************
c calculates Energy and Forces of HLink s
c ramber=rclas
	subroutine link2(numlink,linkat,linkqm,linkmm,linkmm2,
     .    ramber,natot,na_u,namber,fdummy,ng1,attype,nparm,
     .    nbond,nangle,ndihe,nimp,multidihe,multiimp,kbond,bondeq,
     .    kangle,angleeq,kdihe,diheeq,kimp,impeq,perdihe,perimp,
     .    bondtype,angletype,dihetype,imptype,linkqmtype,
     .    bondxat,Ener,parametro,istp)
	
	implicit none
	integer numlink,linkat(15),linkqm(15,4),linkmm(15,4),
     .  na_u,natot,namber,linkmm2(15,4,3),
     .  ng1(namber,6),bondxat(namber)
	integer i,j,j2,k,l,m,jm,jq,at1,at2,at3,at4,istp,p
	double precision ramber(3,natot),flink(3,natot),
     .  fdummy(3,natot),fce(12),Elink(15),scaling
       integer   nbond,nangle,ndihe,nimp,nparm,multidihe(nparm),
     .    multiimp(nparm)
       double precision   kbond(nparm),bondeq(nparm),kangle(nparm),
     .   angleeq(nparm),kdihe(nparm),diheeq(nparm),kimp(nparm),
     .   impeq(nparm),perdihe(nparm),perdihe2(nparm),perimp(nparm)
       character   bondtype(nparm)*5,angletype(nparm)*8,
     .   dihetype(nparm)*11,imptype(nparm)*11,
     .   attype(namber)*4,linkqmtype(15,4)*4
       character ty1*4,ty2*4,ty3*4,ty4*4,tybond*5,
     . tyangle*8,tydihe*11 
       double precision  rij,dist,dx,dy,dz,pi,angle,
     . scal,scalar,dscalar,r12,r32,dr12r32,angulo,
     . dihedro,dihedral,fpar(3),fpp(3),fmod,
     . dversor(3),rhq,rqm,Ener,distl(15)
	logical wriok


	integer ii, ji


c numerical variables to indicate the parameters
c last value eq 0, 1 dihedral, eq 1 more than 1 dihedral
	integer parametro(15,22,4)	

c linkmm2 asigantion: linkmm atoms neighbours
        perdihe2=perdihe
        pi=DACOS(-1.d0)
	wriok=.false.

c change units 
	ramber(1:3,1:natot)=ramber(1:3,1:natot)*0.529177d0

c asignation for E and F
        if(istp.eq.1) then          

 	parametro(1:15,1:22,4)=0

c	bond Cqm -- Cmm parameter 1

	do i=1,numlink
         do k=1,nbond
          tybond=bondtype(k)
          ty1=tybond(1:2)
          ty2=tybond(4:5)
          if(linkqmtype(i,1).eq.ty1.and.attype(linkmm(i,1)).eq.ty2) then
          parametro(i,1,1)=k
          elseif(linkqmtype(i,1).eq.ty2.and.
     .  attype(linkmm(i,1)).eq.ty1) then
          parametro(i,1,1)=k
         endif
         enddo
	enddo

c  H  alyphatihic (1) 340.0 1.090 or aromathic (2) 367.0 1.080

c	do i=1,numlink
c
c        if(linkqmtype(i,1).eq.'CT'.or.linkqmtype(i,1).eq.'CF') then
c	kch(i)=340.0
c        rch(i)=1.09	
c	elseif(linkqmtype(i,1).eq.'CA'.or.linkqmtype(i,1).eq.'CB'.
c     .	or.linkqmtype(i,1).eq.'CC') then
c	kch(i)=367.0
c	rch(i)=1.08
c	else
c	write(*,*) 'Wrong Hlink type  ',linkqmtype(i,1)
c	STOP
c	endif
c	enddo

c 	angles Cqm -- Cmm -- X parameters 2 to 4

	do i=1,numlink
	do j=2,4
	if(linkmm(i,j).eq.0) goto 30
	do k=1,nangle
          tyangle=angletype(k)
          ty1=tyangle(1:2)
          ty2=tyangle(4:5)
          ty3=tyangle(7:8)
          if(linkqmtype(i,1).eq.ty1.and.attype(linkmm(i,1)).eq.ty2.and.
     .       attype(linkmm(i,j)).eq.ty3) then
          parametro(i,j,1)=k
	  elseif(linkqmtype(i,1).eq.ty3.and.attype(linkmm(i,1)).eq.ty2.
     .    and.attype(linkmm(i,j)).eq.ty1) then
          parametro(i,j,1)=k	
	  endif
	enddo
 30	enddo
	enddo
	
c	dihedral X--Cq -- Cmm--Y parameters 5 to 13

	do i=1,numlink
        do jq=2,4
c no more neighbours -> out 
	if(linkqm(i,jq).eq.0) goto 40	
	
	do jm=2,4
c only 2 mm neigh -> out
	if(linkmm(i,jm).eq.0) goto 39

         m=0
         do k=1,ndihe
         tydihe=dihetype(k)
         ty1=tydihe(1:2)
         ty2=tydihe(4:5)
         ty3=tydihe(7:8)
         ty4=tydihe(10:11)
c if X the same constant
         if(ty1.eq.'X ') then
         if(linkqmtype(i,1).eq.ty2.and.
     .      attype(linkmm(i,1)).eq.ty3)  then
         parametro(i,4+3*(jq-2)+jm-1,1)=k
         elseif(linkqmtype(i,1).eq.ty3.and.
     .      attype(linkmm(i,1)).eq.ty2)  then
         parametro(i,4+3*(jq-2)+jm-1,1)=k
         endif
c if not X different constants 
         elseif(ty1.ne.'X ') then

         if(linkqmtype(i,jq).eq.ty1.and.linkqmtype(i,1).eq.
     .   ty2.and.attype(linkmm(i,1)).eq.ty3.and.
     .   attype(linkmm(i,jm)).eq.ty4) then
         parametro(i,4+3*(jq-2)+jm-1,1)=k
   
        if(perdihe2(k).lt.0) then
          m=m+1
	parametro(i,4+3*(jq-2)+jm-1,m)=k           
        parametro(i,4+3*(jq-2)+jm-1,m+1)=k+1   
          endif
         elseif(linkqmtype(i,jq).eq.ty4.and.linkqmtype(i,1).eq.
     .   ty3.and.attype(linkmm(i,1)).eq.ty2.and.
     .   attype(linkmm(i,jm)).eq.ty1) then
	 parametro(i,4+3*(jq-2)+jm-1,1)=k         

          if(perdihe2(k).lt.0) then
          m=m+1
        parametro(i,4+3*(jq-2)+jm-1,m)=k
        parametro(i,4+3*(jq-2)+jm-1,m+1)=k+1 
	 endif
         endif
         endif
         enddo
   39    enddo

   40   enddo
	enddo
	
c	dihedral Cq--Cmm--Y--Z parameters 14 to 22

	do i=1,numlink
        do jm=2,4
	if(linkmm(i,jm).eq.0) goto 51
	do j=1,3
	if(linkmm2(i,jm,j).eq.0) then
	goto 50
	else
         m=0
         do k=1,ndihe
         tydihe=dihetype(k)
         ty1=tydihe(1:2)
         ty2=tydihe(4:5)
         ty3=tydihe(7:8)
         ty4=tydihe(10:11)
c if X the same constant     
         if(ty1.eq.'X ') then
         if(attype(linkmm(i,1)).eq.ty2.and.
     .      attype(linkmm(i,jm)).eq.ty3)  then
         parametro(i,13+3*(jm-2)+j,1)=k
         elseif(attype(linkmm(i,1)).eq.ty3.and.
     .      attype(linkmm(i,jm)).eq.ty2)  then
         parametro(i,13+3*(jm-2)+j,1)=k
         endif
c if not X different constants     
         elseif(ty1.ne.'X ') then
 
         if(linkqmtype(i,1).eq.ty1.and.attype(linkmm(i,1)).eq.
     .   ty2.and.attype(linkmm(i,jm)).eq.ty3.and.
     .   attype(linkmm2(i,jm,j)).eq.ty4) then
	 
        if(perdihe2(k).lt.0) then
        parametro(i,13+3*(jm-2)+j,1)=k
 	 if(perdihe2(k+1).lt.0) then
         parametro(i,13+3*(jm-2)+j,2)=k+1
	  if(perdihe2(k+2).lt.0) then
          parametro(i,13+3*(jm-2)+j,3)=k+2
          parametro(i,13+3*(jm-2)+j,4)=k+3
	  else
	  parametro(i,13+3*(jm-2)+j,3)=k+2
	  goto 50
          endif
         else
	 parametro(i,13+3*(jm-2)+j,2)=k+1
	 goto 50
	 endif
	else	  
	parametro(i,13+3*(jm-2)+j,1)=k
	goto 50
        endif

        elseif(linkqmtype(i,1).eq.ty4.and.attype(linkmm(i,1)).eq.
     .   ty3.and.attype(linkmm(i,jm)).eq.ty2.and.
     .   attype(linkmm2(i,jm,j)).eq.ty1) then 
 
        if(perdihe2(k).lt.0) then
        parametro(i,13+3*(jm-2)+j,1)=k
         if(perdihe2(k+1).lt.0) then
         parametro(i,13+3*(jm-2)+j,2)=k+1
          if(perdihe2(k+2).lt.0) then
          parametro(i,13+3*(jm-2)+j,3)=k+2
           parametro(i,13+3*(jm-2)+j,4)=k+3
          else
          parametro(i,13+3*(jm-2)+j,3)=k+2
          goto 50
          endif
         else
         parametro(i,13+3*(jm-2)+j,2)=k+1
         goto 50
         endif
        else
        parametro(i,13+3*(jm-2)+j,1)=k
        goto 50	      	 
        endif

c one side or other
         endif
c with or without X
         endif
c end dihedrals 
         enddo
	endif
c enddo for neigh 2, neigh 1, numlink 
   50   enddo
   51   enddo
c end numlink
	enddo

c end asignation

c	some printting for debugging	
	if(wriok) then

	write(*,*) 'LINK ATOM VARIABLES'

c neigh matrix

	write(*,*)''
	do i=1,numlink
	write(*,*) 'Neigh Matrix: ',i
	write(*,*) 'Cq:',linkqm(i,1),'Vec: ',(linkqm(i,k),k=2,4)
	write(*,*) 'Cm:	  	    ',linkmm(i,1)
	write(*,*) 'Neigh CM:   ',(linkmm(i,k),k=2,4)
	write(*,*) 'Neigh CM 2: ',((linkmm2(i,j,k),k=1,3),j=2,4)	
	write(*,*)''

	
c cq-cm
	
	write(*,*) 'bond Cq-Cm: ',linkqm(i,1),linkmm(i,1)
	write(*,*) 'types: ',linkqmtype(i,1),attype(linkmm(i,1))
	write(*,*) 'Bondtypenumber(k) ',parametro(i,1,1),
     .   kbond(parametro(i,1,1))
	write(*,*)'************************************'

c Cqm -- Cmm -- X parameters 2 to 4

	write(*,*) 'Angle Cq-Cm-x: ',linkqm(i,1),linkmm(i,1)
	write(*,*)''
	write(*,*) ' x: ',linkmm(i,2),linkmm(i,3),linkmm(i,4)
        write(*,*) 'x types: ',attype(linkmm(i,2)),
     .   attype(linkmm(i,3)),attype(linkmm(i,4))
        write(*,'(a12,i4,f9.3,i4,f9.3,i4,f9.3)')
     .       'Angtype(k): ',
     .	parametro(i,2,1),kangle(parametro(i,2,1)),
     .  parametro(i,3,1),kangle(parametro(i,3,1)),
     .  parametro(i,4,1),kangle(parametro(i,4,1))
	write(*,*)'***********************************'
	
c X--Cq -- Cmm--Y parameters 5 to 13     

	write(*,*)''
	write(*,*) 'Dihedrals y-Cq-Cm-x: ',linkqm(i,1),linkmm(i,1)
	write(*,*)'' 
        write(*,*) ' x: ',linkqm(i,2),linkqm(i,3),linkqm(i,4)
	write(*,*) 'x types: ',linkqmtype(i,2),
     .        linkqmtype(i,3),linkqmtype(i,4)	
	write(*,*)''
	write(*,*) ' y: ',linkmm(i,2),linkmm(i,3),linkmm(i,4)
        write(*,*) 'y types: ',attype(linkmm(i,2)),
     .       attype(linkmm(i,3)),attype(linkmm(i,4))
	write(*,*)''
c x=1
        write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')
     .       'dihe x=1,y=1: (mult) ',
     .       (parametro(i,5,k),kdihe(parametro(i,5,k)),k=1,4)
        write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')
     .       'dihe x=1,y=2: (mult) ',
     .       (parametro(i,6,k),kdihe(parametro(i,6,k)),k=1,4)
        write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')
     .       'dihe x=1,y=3: (mult) ',
     .       (parametro(i,7,k),kdihe(parametro(i,7,k)),k=1,4)
	write(*,*)''
c x=2
        write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')
     .       'dihe x=2,y=1: (mult) ',
     .       (parametro(i,8,k),kdihe(parametro(i,8,k)),k=1,4)
        write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')
     .       'dihe x=2,y=2: (mult) ',
     .       (parametro(i,9,k),kdihe(parametro(i,9,k)),k=1,4)
        write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')
     .       'dihe x=2,y=3: (mult) ',
     .       (parametro(i,10,k),kdihe(parametro(i,10,k)),k=1,4)
	write(*,*)''
c x=3
        write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')
     .       'dihe x=3,y=1: (mult) ',
     .       (parametro(i,11,k),kdihe(parametro(i,11,k)),k=1,4)
        write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')
     .       'dihe x=3,y=2: (mult) ',
     .       (parametro(i,12,k),kdihe(parametro(i,12,k)),k=1,4)
        write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')
     .       'dihe x=3,y=3: (mult) ',
     .       (parametro(i,13,k),kdihe(parametro(i,13,k)),k=1,4)
	write(*,*)'*********************************************'

c Cq--Cmm--Y--Z parameters 14 to 22 

	write(*,*) 'Dihedrals Cq-Cm-y-z: ',linkqm(i,1),linkmm(i,1)
	write(*,*)''
        write(*,*) ' y: ',linkmm(i,2),linkmm(i,3),linkmm(i,4)
        write(*,*) 'y types: ',attype(linkmm(i,2)),
     .        attype(linkmm(i,3)),attype(linkmm(i,4))
 
	write(*,*)''

	do k=1,3
	if(linkmm2(i,2,k).eq.0) then
	goto 60
	else
	l=13+k
        write(*,*) ' z:(y1) ',linkmm2(i,2,k)
        write(*,*) 'z(y1) types: ',attype(linkmm2(i,2,k))
	write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')
     .       'dihedro : (mult)   ',
     .       (parametro(i,l,m),kdihe(parametro(i,l,m)),m=1,4)
	endif
 60	enddo

	do k=1,3
        if(linkmm2(i,3,k).eq.0) then
        goto 61
        else
        l=16+k

	write(*,*) ' z:(y2) ',linkmm2(i,3,k)
        write(*,*) 'z(y2) types: ',attype(linkmm2(i,3,k))
        write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')
     .       'dihedro : (mult)   ',
     .       (parametro(i,l,m),kdihe(parametro(i,l,m)),m=1,4)
	endif
 61	enddo

	do k=1,3
        if(linkmm2(i,4,k).eq.0) then
        goto 62
        else
        l=19+k

        write(*,*) ' z:(y3) ',linkmm2(i,4,k)
        write(*,*) 'z(y3) types: ',attype(linkmm2(i,4,k))
        write(*,'(a20,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3,i2,f6.3)')
     .       'dihedro : (mult)   ',
     .       (parametro(i,l,m),kdihe(parametro(i,l,m)),m=1,4)

	endif
 62	enddo
	write(*,*) '*********************************************'

  	enddo
c wriok
	endif 

c istp=1 only
        endif

c reasignation of perdihe2
        do i=1,ndihe
        if(perdihe2(i).lt.0) then
        perdihe2(i)=-1.d0*perdihe2(i)
        endif
        enddo

	Ener=0.d0

	flink(1:3,1:natot)=0.d0

c calculation of E and F 

	do i=1,numlink
	Elink(i)=0.d0
	
c       bond Cqm -- Cmm parameter 1 
	dx=0.d0
        dy=0.d0
        dz=0.d0

	at1=linkqm(i,1)
	at2=natot-namber+linkmm(i,1)
	k=parametro(i,1,1)
	
         rij=dist(ramber(1,at1),ramber(2,at1),ramber(3,at1),
     .   ramber(1,at2),ramber(2,at2),ramber(3,at2))

          Elink(i) = Elink(i) + kbond(k)
     .    *(rij-bondeq(k))**2 
 	
          dx=(1.d0/rij)*(ramber(1,at1)-ramber(1,at2))
          dx=2.d0*kbond(k)*(rij-bondeq(k))*dx
 
          dy=(1.d0/rij)*(ramber(2,at1)-ramber(2,at2))
          dy=2.d0*kbond(k)*(rij-bondeq(k))*dy
 
          dz=(1.d0/rij)*(ramber(3,at1)-ramber(3,at2))
          dz=2.d0*kbond(k)*(rij-bondeq(k))*dz
	 
c QM atoms are < 0, atom 2 are equal 
        flink(1,at1)=-dx
        flink(2,at1)=-dy
        flink(3,at1)=-dz
        flink(1,at2)=dx
        flink(2,at2)=dy
        flink(3,at2)=dz

c        write(*,*) "Hlink bonds"
c        write(*,*) i, Elink(i), flink


c	write(*,*) 'contrib cq-cm i E F',i,
c     .  kbond(k)*(1.0-kbond(k)/kch(i))
c     .    *(rij-bondeq(k))**2,
c     .   sqrt(flink(1,at1)**2+flink(2,at1)**2+flink(3,at1)**2)

c ***********************************************************************
c       angle  Cqm -- Cmm -- X parameters 2 to 4

c	3 angles ( one for each x) j, x index
	dx=0.d0
	dy=0.d0
	dz=0.d0
	
	do j=2,4
	if(linkmm(i,j).eq.0) goto 55
	at1=linkqm(i,1)
	at2=natot-namber+linkmm(i,1)
	at3=natot-namber+linkmm(i,j)
	k=parametro(i,j,1)

	        angulo=angle(ramber(1,at1),ramber(2,at1),ramber(3,at1),
     .          ramber(1,at2),ramber(2,at2),ramber(3,at2),
     .          ramber(1,at3),ramber(2,at3),ramber(3,at3))
 
         Elink(i) = Elink(i) + kangle(k)*
     .                ((angulo-angleeq(k))*
     .                (pi/180d0))**2d0
 
c	write(*,*) 'contrib cq-cm-x(j) i,j E ',i,j,kangle(k)*
c     .                ((angulo-angleeq(k))*
c     .                (pi/180))**2
 	
        scal=scalar(ramber(1,at1),ramber(2,at1),ramber(3,at1),
     .          ramber(1,at2),ramber(2,at2),ramber(3,at2),
     .          ramber(1,at3),ramber(2,at3),ramber(3,at3))
 
        r12=dist(ramber(1,at1),ramber(2,at1),ramber(3,at1),
     .  ramber(1,at2),ramber(2,at2),ramber(3,at2))
 
        r32=dist(ramber(1,at3),ramber(2,at3),ramber(3,at3),
     .  ramber(1,at2),ramber(2,at2),ramber(3,at2))
	      
	dscalar=(ramber(1,at3)-ramber(1,at2))
        dr12r32=r32*(ramber(1,at1)-ramber(1,at2))/(r12)
        dx=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
        dx=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dx
        dx=2.d0*kangle(k)*(angulo-angleeq(k))*
     .                (pi/180d0)*dx
 
        dscalar=(ramber(2,at3)-ramber(2,at2))
        dr12r32=r32*(ramber(2,at1)-ramber(2,at2))/(r12)
        dy=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
        dy=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dy
        dy=2.d0*kangle(k)*(angulo-angleeq(k))*
     .                (pi/180d0)*dy
 
        dscalar=(ramber(3,at3)-ramber(3,at2))
        dr12r32=r32*(ramber(3,at1)-ramber(3,at2))/(r12)
        dz=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
        dz=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dz
        dz=2.d0*kangle(k)*(angulo-angleeq(k))*
     .                (pi/180d0)*dz
	
c atom qm cq force, sum of each particular x
	flink(1,at1)=flink(1,at1)-dx
	flink(2,at1)=flink(2,at1)-dy
	flink(3,at1)=flink(3,at1)-dz
	
c force over x, change at3 by at1 
	dx=0.d0
        dy=0.d0
        dz=0.d0
	dscalar=(ramber(1,at1)-ramber(1,at2))
        dr12r32=r32*(ramber(1,at3)-ramber(1,at2))/(r12)
        dx=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
        dx=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dx
        dx=2.d0*kangle(k)*(angulo-angleeq(k))*
     .                (pi/180d0)*dx
 
        dscalar=(ramber(2,at1)-ramber(2,at2))
        dr12r32=r32*(ramber(2,at3)-ramber(2,at2))/(r12)
        dy=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
        dy=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dy
        dy=2.d0*kangle(k)*(angulo-angleeq(k))*
     .                (pi/180d0)*dy
 
        dscalar=(ramber(3,at1)-ramber(3,at2))
        dr12r32=r32*(ramber(3,at3)-ramber(3,at2))/(r12)
        dz=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
        dz=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dz
        dz=2.d0*kangle(k)*(angulo-angleeq(k))*
     .                (pi/180d0)*dz
	

c sum each x  
	flink(1,at3)=flink(1,at3)-dx
        flink(2,at3)=flink(2,at3)-dy
        flink(3,at3)=flink(3,at3)-dz
	
c middle atom force (change in derivate) 
	dx=0.d0
        dy=0.d0
        dz=0.d0
	
      dscalar=2.d0*ramber(1,at2)-ramber(1,at1)-
     .        ramber(1,at3)
      dr12r32=(r32*(-ramber(1,at1)+ramber(1,at2))/r12)+
     .        (r12*(-ramber(1,at3)+ramber(1,at2))/r32)
      dx=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
      dx=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dx
      dx=2.d0*kangle(k)*(angulo-angleeq(k))*
     .                (pi/180d0)*dx
 
      dscalar=2.d0*ramber(2,at2)-ramber(2,at1)-
     .        ramber(2,at3)
      dr12r32=(r32*(-ramber(2,at1)+ramber(2,at2))/r12)+
     .        (r12*(-ramber(2,at3)+ramber(2,at2))/r32)
      dy=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
      dy=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dy
      dy=2.d0*kangle(k)*(angulo-angleeq(k))*
     .                (pi/180d0)*dy
 
      dscalar=2.d0*ramber(3,at2)-ramber(3,at1)-
     .        ramber(3,at3)
      dr12r32=(r32*(-ramber(3,at1)+ramber(3,at2))/r12)+
     .        (r12*(-ramber(3,at3)+ramber(3,at2))/r32)
      dz=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
      dz=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dz
      dz=2.d0*kangle(k)*(angulo-angleeq(k))*
     .                (pi/180d0)*dz

c central atom force 
        flink(1,at2)=flink(1,at2)-dx
        flink(2,at2)=flink(2,at2)-dy
        flink(3,at2)=flink(3,at2)-dz

c enddo for each possible x 
 55	enddo

c        write(*,*) "Hlink angulos"
c        write(*,*) i, Elink(i), flink



c********************************************************
c	dihedrals  X--Cq -- Cmm--Y parameters 5 to 13
c	for each 3 X 
c	for each 3 Y
	do jq=2,4
	do jm=2,4
	
	if(linkmm(i,jm).eq.0) goto 71
	at1=linkqm(i,jq)
	at2=linkqm(i,1)   
	at3=natot-namber+linkmm(i,1)
	at4=natot-namber+linkmm(i,jm)      

c	write(*,*) "flag Js", jq,jm,at1, at2, at3,at4

	dihedral=0
	do m=1,4
	k=parametro(i,4+3*(jq-2)+jm-1,m)
	  if (k.ne.0) then 
c modificacion para debug, Nick
c antes queria calcular diedros usando un atomo 0
c y luego mataba el calculo con k=0
	dihedral=dihedro(ramber(1,at1),ramber(2,at1),ramber(3,at1),
     .   ramber(1,at2),ramber(2,at2),
     .   ramber(3,at2),
     .   ramber(1,at3),ramber(2,at3),
     .   ramber(3,at3),
     .   ramber(1,at4),ramber(2,at4),
     .   ramber(3,at4))
c	  else
c	    dihedral=0
	  end if
	end do

c	write(*,*) "diedro", dihedral

c who many parameters for dihedrals 5 to 13
	do m=1,4		
        k=parametro(i,4+3*(jq-2)+jm-1,m)
	if(k.ne.0) then	

c	write(*,*) "Emod", k, kdihe(k), multidihe(k), perdihe2(k),diheeq(k)

      Elink(i) =Elink(i)+(kdihe(k)/multidihe(k))*
     .  (1+dCOS((pi/180d0)*(perdihe2(k)*dihedral-diheeq(k))))
c	write(*,*) "Elink", i,jq,jm,m,Elink(i)
        
c	write(*,*) 'contrib x-cq-cm-y ',i,jq,jm,(kdihe(k)/multidihe(k))*
c     .  (1+dCOS((pi/180)*(perdihe2(k)*dihedral-diheeq(k))))
	
c force for each 4 atoms 
	do j=1,4
        call diheforce(natot,ramber,
!cmabie namber por natot, parece q sino esta mal la dimension de ramber en diheforce, Nick
     .                 at1,at2,at3,at4,j,
     .                 kdihe(k),diheeq(k),perdihe2(k),multidihe(k),fce)

c force assignation 
	flink(1,at1)=flink(1,at1)-fce(1)
        flink(2,at1)=flink(2,at1)-fce(2)
        flink(3,at1)=flink(3,at1)-fce(3) 
	flink(1,at2)=flink(1,at2)-fce(4)
        flink(2,at2)=flink(2,at2)-fce(5)
        flink(3,at2)=flink(3,at2)-fce(6)
        flink(1,at3)=flink(1,at3)-fce(7)
        flink(2,at3)=flink(2,at3)-fce(8)
        flink(3,at3)=flink(3,at3)-fce(9)
	flink(1,at4)=flink(1,at4)-fce(10)
        flink(2,at4)=flink(2,at4)-fce(11)
        flink(3,at4)=flink(3,at4)-fce(12)	

	enddo
	else
	goto 70
	endif
c enddo more than 1 parameter 
   70	enddo

c enddo X and Y 
  71	enddo
	enddo
	

c 	dihedrals Cq -- Cmm--Y--Z  parameters 14 to 22
c 	for each 3 X
c 	for each 3 Y
        
        do jm=2,4
	if(linkmm(i,jm).eq.0) goto 81
	do jq=1,3 
c is other Z atom? 
	if(linkmm2(i,jm,jq).eq.0) then
        goto 80
	else
	l=13+3*(jm-2)+jq
        at1=linkqm(i,1)
        at2=natot-namber+linkmm(i,1)
        at3=natot-namber+linkmm(i,jm)
	at4=natot-namber+linkmm2(i,jm,jq) 
 
        dihedral=dihedro(ramber(1,at1),ramber(2,at1),ramber(3,at1),
     .   ramber(1,at2),ramber(2,at2),
     .   ramber(3,at2),
     .   ramber(1,at3),ramber(2,at3),
     .   ramber(3,at3),
     .   ramber(1,at4),ramber(2,at4),
     .   ramber(3,at4))


c how many parameters for dihedrals 14 to 22 
c more than 1 paramter
        do m=1,4
        k=parametro(i,l,m)
        if(k.ne.0) then
 
      Elink(i) =Elink(i)+(kdihe(k)/multidihe(k))*
     .  (1+dCOS((pi/180d0)*(perdihe2(k)*dihedral-diheeq(k))))
 
c force over 4 atoms
        do j=1,4
        call diheforce(natot,ramber,
!idem antes namber -> natot
     .                 at1,at2,at3,at4,j,
     .                 kdihe(k),diheeq(k),perdihe2(k),multidihe(k),fce)


c force assignation 
        flink(1,at1)=flink(1,at1)-fce(1)
        flink(2,at1)=flink(2,at1)-fce(2)
        flink(3,at1)=flink(3,at1)-fce(3)
        flink(1,at2)=flink(1,at2)-fce(4)
        flink(2,at2)=flink(2,at2)-fce(5)
        flink(3,at2)=flink(3,at2)-fce(6)
        flink(1,at3)=flink(1,at3)-fce(7)
        flink(2,at3)=flink(2,at3)-fce(8)
        flink(3,at3)=flink(3,at3)-fce(9)
        flink(1,at4)=flink(1,at4)-fce(10)
        flink(2,at4)=flink(2,at4)-fce(11)
        flink(3,at4)=flink(3,at4)-fce(12)
	enddo
	
        else
	goto 75
	endif
c enddo more than 1 parameter 
   75   enddo
	 
	endif
c enddo to Y and Z
   80   enddo
   81   enddo


c ***************************************************
c For each LA Force over HL scaling and sum where corresponds 
c HL force division in parallel fpar and perpendicular fpp
c unitary versor dversor(3) at3=LA
	scaling=0.1
	at1=linkqm(i,1)
	at2=natot-namber+linkmm(i,1)
	at3=linkat(i)
	
	dversor(1)=ramber(1,at1)-ramber(1,at3)
	dversor(2)=ramber(2,at1)-ramber(2,at3)
	dversor(3)=ramber(3,at1)-ramber(3,at3)
	
         rij=dist(ramber(1,at1),ramber(2,at1),ramber(3,at1),
     .   ramber(1,at3),ramber(2,at3),ramber(3,at3))
	
	dversor(1:3)=dversor(1:3)/rij

c angle proyection calculation, to proyect forces 
	angulo=0.d0
       angulo=angle(fdummy(1,at3),fdummy(2,at3),fdummy(3,at3),
     .          0.d0,0.d0,0.d0,
     .          dversor(1),dversor(2),dversor(3) )

	angulo=angulo*pi/180.d0
c        write(*,*) "angulo", angulo

c  product  modF and dversor
	 fmod=sqrt(fdummy(1,at3)**2d0+
     .   fdummy(2,at3)**2d0+fdummy(3,at3)**2d0)

	fpar(1:3)=fmod*dversor(1:3)*dCOS(angulo)
c	write(*,*) "fpar calc", fmod, dversor(1:3), dCOS(angulo)
	fpp(1:3)=fdummy(1:3,at3)-fpar(1:3) 
c	write(*,*) "fpp calc", fdummy(1:3,at3),fpar(1:3)

c perpendicular fce division and sum to CMM y CQM
c parallel fce sum over CM and HL fce eq cero
	fdummy(1:3,at3)=fpar(1:3)*scaling

c distance HL--CQ and CM--CQ
	rhq=dist(ramber(1,at1),ramber(2,at1),ramber(3,at1),
     .  ramber(1,at3),ramber(2,at3),ramber(3,at3))

	rqm=dist(ramber(1,at1),ramber(2,at1),ramber(3,at1),
     .  ramber(1,at2),ramber(2,at2),ramber(3,at2))
	


	flink(1,at2)=flink(1,at2)+fpp(1)*rhq/rqm
        flink(2,at2)=flink(2,at2)+fpp(2)*rhq/rqm 
        flink(3,at2)=flink(3,at2)+fpp(3)*rhq/rqm 

	flink(1,at1)=flink(1,at1)+fpp(1)*(1.d0-rhq/rqm)
        flink(2,at1)=flink(2,at1)+fpp(2)*(1.d0-rhq/rqm)
        flink(3,at1)=flink(3,at1)+fpp(3)*(1.d0-rhq/rqm)


	Ener=Ener+Elink(i)

c      enddo numlink
 	enddo	


	fdummy(1:3,1:natot)=fdummy(1:3,1:natot)+flink(1:3,1:natot)    

c	write(*,*) 'Total LA energy: ',Ener 

c change units 
        ramber(1:3,1:natot)=ramber(1:3,1:natot)/0.529177d0

	end

c*****************************************************************
c calculates HL position 
	
	subroutine link3(numlink,linkat,linkqm,linkmm,ramber,
     .    natot,na_u,namber,distl)

        implicit none
        integer numlink,linkat(15),linkqm(15,4),linkmm(15,4),
     .  na_u,natot,namber,linkmm2(15,4,3)
     
        integer i,j,j2,k,l,m,jm,jq,at1,at2,at3,at4
        double precision ramber(3,natot),
     .	x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
     .  r4,a4,d4,angle,dihedro2,distl(15)
        logical frstme
        save frstme
        data frstme /.true./

c change units
        ramber(1:3,1:natot)=ramber(1:3,1:natot)*0.529177d0

C dist HL-CQM  
	do i=1,numlink
       distl(i)=sqrt((ramber(1,linkqm(i,1))-ramber(1,linkat(i)))**2d0
     . +(ramber(2,linkqm(i,1))-ramber(2,linkat(i)))**2d0
     . +(ramber(3,linkqm(i,1))-ramber(3,linkat(i)))**2d0)
	enddo

C sets HL-CQM dist first time
	if(frstme) then
        frstme=.false.
	distl=1.09
	endif

C sets HL new positions
	do i=1,numlink
	at2=linkqm(i,1)
	at1=natot-namber+linkmm(i,1)
	at3=linkqm(i,2)
	at4=linkqm(i,3)

c angle CQM-CMM
	a4=angle(ramber(1,at1),ramber(2,at1),ramber(3,at1),
     .          ramber(1,at2),ramber(2,at2),ramber(3,at2),
     .          ramber(1,at3),ramber(2,at3),ramber(3,at3))

c diehdro CQM1-HL-CQM2-CQM3
        d4=dihedro2(ramber(1,at1),ramber(2,at1),ramber(3,at1),
     .   ramber(1,at2),ramber(2,at2),
     .   ramber(3,at2),
     .   ramber(1,at3),ramber(2,at3),
     .   ramber(3,at3),
     .   ramber(1,at4),ramber(2,at4),
     .   ramber(3,at4))

c HL-CQM distance 
	r4=distl(i)

	x1=ramber(1,at4)
	y1=ramber(2,at4)
	z1=ramber(3,at4)
        x2=ramber(1,at3)
        y2=ramber(2,at3)
        z2=ramber(3,at3)
        x3=ramber(1,at2)
        y3=ramber(2,at2)
        z3=ramber(3,at2)

c change HL position
	call pos4(x1,y1,z1,x2,y2,z2,x3,y3,z3,r4,a4,d4,x4,y4,z4)  

c HL coords change 
	ramber(1,linkat(i))=x4
	ramber(2,linkat(i))=y4
	ramber(3,linkat(i))=z4

c enndo numlink	
	enddo

c change units
        ramber(1:3,1:natot)=ramber(1:3,1:natot)/0.529177d0

	end

c *************************************************************************
