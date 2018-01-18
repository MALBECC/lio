c subroutine that read the constrained atom block 

	subroutine fixed1(na_u,nac,natot,nroaa,rclas,blocklist,
     .                    atname,aaname,aanum,wat)

        use precision
        use fdf
        use sys
	implicit none
	integer i,j,k,l,iunit,na_u,nac,natot,nroaa,aanum(nac),
     .  blocklist(natot)
	character*4 atname(nac),aaname(nac)
	double precision rclas(3,natot)
	logical wat

	integer nt,type(9),aa,con1(4,natot),con2(4,natot),ncon(4)
	integer nc,c1(na_u),c2(na_u)
	double precision r(3,nac),cut,cqm(3),mdist,dist,dist2
	character ch1*1,ch4*4,exp
	wat = .false.

c reads 'PositionConstraints' block
	if ( fdf_block('PositionConstraints',iunit) ) then
        read(iunit,'(A)',advance='no',err=100,end=100) exp 
        if(exp.eq.'%') goto 50 
	read(iunit,*,err=100,end=100) exp,nt
	do i=1,nt
	  read(iunit,*,err=100,end=100) exp,type(i)

	if(nac.eq.0.and.type(i).gt.1) then
	call die('fixed: constrained type for only QM atoms
     .  must be only 1')
	endif
	  if(type(i).le.2) then

       k=1
 10    continue
       if(k.gt.natot) then
       call die('fixed: sets of constraints must not exeed natot')
       endif
       read(iunit,'(A)',advance='no',err=100,end=100) exp
       if(exp.eq.'f'.or.exp.eq.'F'.or.exp.eq.'p'.or.exp.eq.'P') then
       read(iunit,*,err=100,end=100) exp,con1(i,k),exp,con2(i,k)
       k=k+1
       goto 10 
       else
       continue
       endif
       ncon(i)=k-1

	if(ncon(i).eq.0) then
	call die('fixed: sets of constraints must not be zero') 
	endif

	do k=1,ncon(i)
	if(con2(i,k).lt.con1(i,k)) then
	call die('fixed: sets of constraints must be defined
     .  in correct order')
	endif
	enddo

	  elseif(type(i).eq.6) then
            read(iunit,*,err=100,end=100) exp,aa
            read(iunit,*,err=100,end=100) exp,cut
	  elseif(type(i).lt.-2 .or. type(i).gt.7) then
	  call die('fixed: Wrong type of constraint')
	  endif !type
	enddo !nt
	write(6,'(/,a)') 'fixed: reading "PositionConstraints" block'
	else 
	goto 50
	endif !fdf

c blocklist assignation
      do i=1,nt
	if(type(i).eq.1) then
          do k=1,ncon(i)
        if(con2(i,k).gt.natot) then
        call die('fixed: atoms in constraint must not exeed natot')
        endif
            do j=con1(i,k),con2(i,k)
            blocklist(j)=1
            enddo
          enddo
	elseif(type(i).eq.-1) then
        blocklist=1
          do k=1,ncon(i)
        if(con2(i,k).gt.natot) then
        call die('fixed: atoms in constraint must not exeed natot')
        endif
            do j=con1(i,k),con2(i,k)
            blocklist(j)=0
            enddo
          enddo
	elseif(type(i).eq.2) then
          do k=1,ncon(i)
        if(con2(i,k).gt.nroaa) then
        call die('fixed: residues in constraint must not exeed nroaa')
        endif

            do j=con1(i,k),con2(i,k)
              do l=1,nac
                if(j.eq.aanum(l)) then
                blocklist(na_u+l)=1
                endif
              enddo
            enddo
          enddo
	elseif(type(i).eq.-2) then
        blocklist=1
          do k=1,ncon(i)
        if(con2(i,k).gt.nroaa) then
        call die('fixed: residues in constraint must not exeed nroaa')
        endif
            do j=con1(i,k),con2(i,k)
              do l=1,nac
                if(j.eq.aanum(l)) then
                blocklist(na_u+l)=0
                endif
              enddo
            enddo
          enddo
	elseif(type(i).eq.3) then
          do l=1,nac
            ch4=atname(l)
            ch1=ch4(1:1)
            if(ch1.ne.'H') then
            blocklist(na_u+l)=1
            endif
          enddo
	elseif(type(i).eq.4) then
          do l=1,nac
            if(atname(l).eq.'CA'.or.
     .         atname(l).eq.'C' .or.
     .         atname(l).eq.'N') then
            blocklist(na_u+l)=1
            endif
          enddo
	elseif(type(i).eq.5) then
          do l=1,nac
            if(aaname(l).ne.'HOH') then
            blocklist(na_u+l)=1
            endif
          enddo
	elseif(type(i).eq.6) then
        r(1:3,1:nac)=rclas(1:3,na_u+1:natot)*0.529177d0
        cqm=0.0
        k=0
        do j=1,nac
        if(aanum(j).eq.aa) then
        k=k+1
        cqm(1)=cqm(1)+r(1,j)
        cqm(2)=cqm(2)+r(2,j)
        cqm(3)=cqm(3)+r(3,j)
        endif
        enddo
        cqm(1:3)=cqm(1:3)/k

        mdist=0.0
        dist=0.0
        do j=1,nac
        if(aanum(j).eq.aa) then
        dist=(r(1,j)-cqm(1))**2+
     .       (r(2,j)-cqm(2))**2+
     .       (r(3,j)-cqm(3))**2
        if(dist.gt.mdist) mdist=dist
        endif
        enddo
        mdist=sqrt(mdist)
        cut=cut+mdist

        dist=0.0
        dist2=cut**2
        do j=1,nac
        if(aaname(j).eq.'HOH'.and.atname(j).eq.'O') then
                        dist=(r(1,j)-cqm(1))**2+
     .                       (r(2,j)-cqm(2))**2+
     .                       (r(3,j)-cqm(3))**2
                        if(dist.gt.dist2) then
                        blocklist(j+na_u)=1
                        blocklist(j+1+na_u)=1
                        blocklist(j+2+na_u)=1
                        endif
        elseif(aaname(j).ne.'HOH') then
                        dist=(r(1,j)-cqm(1))**2+
     .                       (r(2,j)-cqm(2))**2+
     .                       (r(3,j)-cqm(3))**2
                        if(dist.gt.dist2) then
                        blocklist(j+na_u)=1
                        endif
        endif
        enddo
	elseif(type(i).eq.7) then
	wat = .true.
	write(6,'(/,a)') 'fixed: Running restraining water cap'
	endif !type
      enddo !nt

 50    continue 
c read 'GeometryConstraints' block
	if ( fdf_block('GeometryConstraints',iunit) ) then
        k=1
 60     continue
        if(k.gt.na_u) then
        call die('fixed: sets of QM constraints must not exeed na_u')
        endif
        read(iunit,'(A)',advance='no',err=200,end=200) exp
        if(exp.eq.'%'.and.k.eq.1) goto 70
        if(exp.eq.'%'.and.k.ne.1) goto 65
        if(exp.eq.'P'.or.exp.eq.'p') then
        read(iunit,*,err=200,end=200) exp,exp,c1(k),exp,c2(k)
        k=k+1
        goto 60
        else
        call die('fixed: wrong syntaxis in GeometryConstraints block')
        endif
 65     continue
        nc=k-1
 
        do i=1,nc
        if(c2(i).gt.na_u) then
        call die('fixed: atoms in QM constraint must not exeed na_u')
        endif

        if(c2(i).lt.c1(i)) then
        call die('fixed: sets of QM constraints must be defined 
     .            in correct order')
        endif

c blocklist assignation
        do j=c1(i),c2(i)
         blocklist(j)=1
        enddo
        enddo

        write(6,'(/,a)') 'fixed: reading "GeometryConstraints" block'
	endif !fdf 

 70    continue
       return
 100   stop 'fixed: Problem reading "PositionConstraints" block'
 200   stop 'fixed: Problem reading "GeometryConstraints" block'
       end

c*************************************************************
c subroutine that imposes fce and vel constraints 

       subroutine  fixed2(na_u,nac,natot,nfree,blocklist,blockqmmm,
     .             fdummy,cfdummy,vat)

	implicit none
	integer i,k,na_u,nac,natot,nfree,blockqmmm(nac),blocklist(natot)
	double precision fdummy(3,natot),cfdummy(3,natot),vat(3,natot) 
        logical    frstme
        save       frstme
        data frstme /.true./

C copy fdummy array to cfdummy
	cfdummy = fdummy 

c nullify fce and vel 
        do i=1,na_u
           if (blocklist(i).eq.1) then
           cfdummy(1:3,i)=0.0
           vat(1:3,i)=0.0
           endif
        enddo

        do i=1,nac
           if (blockqmmm(i).eq.1.or.blocklist(i+na_u).eq.1) then
           cfdummy(1:3,i+na_u)=0.0
           vat(1:3,i+na_u)=0.0
           endif
        enddo

c set total number of free atoms
	if(frstme) then
	  k=0
	  do i=1,na_u
	    if((blocklist(i).eq.0)) k=k+1
	  enddo
	  do i=1,nac
 	    if((blocklist(i+na_u).eq.0).and.(blockqmmm(i).eq.0)) k=k+1
	  enddo
	  nfree=k
	  write(6,'(/a,2x,i5)') 'hybrid: Total Free Atoms:', nfree 
	  frstme=.false.
	endif

	end
c******************************************************************

