c this subroutine do everything related to reading and parameters
c assignation of solvent atoms

      subroutine solv_assign(na_u,natot,nac,nroaa,Em,Rm,attype,pc,
     .  ng1,bondxat,angexat,atange,angmxat,atangm,dihexat,atdihe,
     .  dihmxat,atdihm,impxat,atimp,
     .  nbond,kbond,bondeq,bondtype,
     .  nangle,kangle,angleeq,angletype,
     .  ndihe,kdihe,diheeq,dihetype,multidihe,perdihe,
     .  nimp,kimp,impeq,imptype,multiimp,perimp,
     .  nparm,aaname,atname,aanum,qmattype,rclas,
     .  rcorteqmmm,rcorteqm,rcortemm,sfc,
     .  radbloqmmm,atsinres, radblommbond)

	use precision
	use sys
	use fdf
	
        implicit none
        integer i,j,k,l,m,n,natot,nac,na_u,nroaa,iunit,                                                              
     .  ncon,nparm,nbond,nangle,ndihe,nimp,atsinres(20000)     
        integer ng1(nac,6),con2(2,1000)
        double precision, dimension(:,:), allocatable, save::
     .  qaa,Rma,Ema 
        double precision Rm(natot),Em(natot),
     .  pc(0:nac),rclas(3,natot)
cagregue 0 apc
        character ch*1,exp
 	character*4 atom
	character*4, dimension(:), allocatable, save::
     .  resname
        character*4, dimension(:,:), allocatable, save::
     .  atnamea,attypea,atsp
	integer, dimension(:,:), allocatable, save::
     .	atnu
        integer, dimension(:,:), allocatable, save::      
     .  nataa,con
	integer, dimension(:), allocatable, save::
     .  atxres
        integer atange(nac,25,2),atangm(nac,25,2),
     .  atdihe(nac,100,3),atdihm(nac,100,3)
        integer bondxat(nac),angexat(nac),
     .  dihexat(nac),dihmxat(nac),angmxat(nac)
        integer  impxat(nac),atimp(nac,25,4)
        character*4 aanamea(200)
        integer atomsxaa(200)
        integer aanum(nac),resnum(nac)
        double precision rcorteqm,rcortemm,sfc,rcorteqmmm
        character*4 atname(nac),aaname(nac),attype(nac),qmattype(na_u)
        character*5 bondtype(nparm)
        character*8 angletype(nparm)
        character*11 dihetype(nparm),imptype(nparm)
        double precision kbond(nparm),bondeq(nparm),kangle(nparm),
     .  angleeq(nparm),kdihe(nparm),diheeq(nparm),perdihe(nparm),
     .  kimp(nparm),impeq(nparm),perimp(nparm)
        integer multidihe(nparm),multiimp(nparm)
        double precision radbloqmmm
	logical foundamber
        character ch1*1,ch4*4
ccorrecion del bug de bonds extras, Nick
	integer :: ivalue(natot*2), inick,jnick
c parche para q no ponga bonds entre extremos terminales
        double precision radblommbond

	ivalue=0


C nullify some solvent vbles
      rclas = 0.d0
      aanum = 0
      nroaa = 0
      ng1 = 0
      nbond = 0
      nangle = 0
      ndihe = 0
      nimp = 0
      ncon = 0
      rcortemm = 100.d0
      rcorteqm = 1.d-06
      rcorteqmmm = 100.d0
      sfc=2.d0
      radbloqmmm=100.d0
      foundamber=.false.

C read solvent coordinates by atom
      if ( fdf_block('SolventInput',iunit) ) then
	do i=1,nac
         read(iunit,err=10,end=10,fmt='(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)')
     .      atom, j, atname(i), aaname(i),
     .      ch, resnum(i), rclas(1:3,na_u+i)

	 ivalue(j)=i
         enddo
	else
        call die("solvent: You must specify the solvent coordinates")      
      endif


c change coordinates to Siesta format
	rclas(1:3,1:natot) = rclas(1:3,1:natot) / 0.529177d0   

c assigns number of residues 
       k=1
       aanum(1)=k
       do i=2,nac
         if (resnum(i).eq.resnum(i-1)) then
           aanum(i)=aanum(i-1)
         elseif (resnum(i).ne.resnum(i-1)) then
           k = k+1
           aanum(i)= k
         endif
       enddo

       nroaa = aanum(nac)

       if(nroaa.eq.0) then
         call die("solvent: Number of residues can not be zero")
       endif
 
       if(nac.eq.0) nroaa=0

c read solute atom type
        if(na_u.ne.0) then
      if ( fdf_block('SoluteAtomTypes',iunit) ) then
	do i=1,na_u+1
	  read(iunit,*,end=20,err=20) ch4 
	  ch1=ch4(1:1)
c	  write(*,*) "ch4 vale ", ch1, " Nick"
	  if(i.eq.na_u+1) then
	    if(ch1.eq.'%') then
	    goto 2
	    else 
	    call die('solvent: solute atom types are greater than na_u') 
	    endif
	  endif
	  if(ch1.eq.'%') then
	  call die('solvent: solute atom types are lower than na_u') 
	  endif
	qmattype(i)=ch4
 2	enddo
      else
       call die('solvent: You must specify solute atom types')
      endif
	endif

c read cutoff radious
      if ( fdf_block('CutOffRadius',iunit) )
     .     then
      read(iunit,*,err=30,end=30) exp, rcorteqm
      read(iunit,*,err=30,end=30) exp, rcorteqmmm
      read(iunit,*,err=30,end=30) exp, rcortemm
      read(iunit,*,err=30,end=30) exp, radbloqmmm
      read(iunit,*,err=30,end=30) exp, radblommbond
C      read(iunit,*,err=30,end=30) exp, sfc
      else
      write(6,'(/a)') 'solvent: Cut-off radius will be the standard'
      endif

c checking cut-off radius
      if(rcorteqm.le.1.e-8) then
      call die('solvent: QM cut-off radius to close to zero')
      endif
      if(rcorteqmmm.le.1.e-8) then
      call die('solvent: QM-MM cut-off radius to close to zero')
      endif

c external solvent connectivity block
       if ( fdf_block('SolventConnectivity',iunit) ) then
         i=1
 5       continue
         if(i.gt.1000) then
           call die('read: Solvent connectivities must not exeed 1000')
         endif
         read(iunit,'(a)',advance='no',err=40,end=40) exp 
         if(exp.eq.'%') goto 6
         read(iunit,*,err=40,end=40) exp, con2(1:2,i)
         i=i+1
         goto 5
 6       continue
         ncon=i-1
         allocate(con(2,ncon))
c arreglo esto para conectividades correctas segun el numero de atomo en el fdf, Nick
c         con(1:2,1:ncon)=con2(1:2,1:ncon)

	  do jnick=1, ncon
	  do inick=1, 2
	     write(*,*) "i,j,con2, ivalue", inick,jnick,con2(inick,jnick),
     .       ivalue(con2(inick,jnick))
	     con(inick,jnick)=ivalue(con2(inick,jnick))
	  end do
	     write(*,*) "interaccion extra agregada entre ",
     .       con(1,jnick), " y ", con(2,jnick)
	  end do

          if(ncon.ne.0) then
            write(6,'(/,a)') 'read: Reading new connectivities block'
          endif
       else
         allocate(con(2,1))
cagregado nick, sino trae problemas cuando ncon=0
       endif


c se fija si esta el amber.parm 
	inquire( file="amber.parm", exist=foundamber )
	if(.not.foundamber) then
	call die("solvent: 'amber.parm' file not found")
	endif

c llama a sub q lee los parametros de bond, angle, dihe e imp
       call amber_union_parms(nbond,kbond,bondeq,bondtype,
     .        nangle,kangle,angleeq,angletype,ndihe,kdihe,
     .        diheeq,dihetype,multidihe,perdihe,nimp,kimp,
     .        impeq,imptype,multiimp,perimp,nparm)

c aca empieza la verdadera asignacion: segun atomo xa c/aa
	allocate(atnu(nroaa,100),atnamea(nroaa,100),
     .	resname(nroaa),atxres(nroaa))

c subrutina que lee el nro de atomos por aa
        call atxaa(n,aanamea,atomsxaa)

c cambia vbles
	atxres=0
        k = 1
        do i=1,nroaa
	   resname(i)=aaname(k)
		do l=1,n
        		if (resname(i).eq.aanamea(l)) then
        		atxres(i)=atomsxaa(l)
        		endif
        	enddo
		if(atxres(i).eq.0) then
			write(6,*) 'solvent: Wrong residue name  :', i
			STOP
		endif

        	do j=1,atxres(i)
        		atnu(i,j)=k
			atnamea(i,j)=atname(k)
        		k = k+1
        	enddo
        enddo	

c llama a la subrutina q trae los parametros amber: attype,pc
	allocate(qaa(nroaa,100),attypea(nroaa,100),nataa(nroaa,100)) 
	call paramats(nroaa,resname,atnamea,qaa,attypea,nataa,
     .                nac,atxres,pc,attype)

c llama a la subrutina q trae los parametros amber: Rm,Em 
	allocate(Ema(nroaa,100),Rma(nroaa,100))
	call lj(nroaa,attypea,Ema,Rma,na_u,natot,Em,Rm,qmattype,atxres) 

c llama a la subrutina q asigna los 1eros vecinos
       call amber_ng1(nac,nataa,nroaa,atxres,atnu,resname,
     .                ng1,atnamea,ncon,con)
 
c llama a la sub q calcula bonds, angulos, dihedros e impropers
       call bon_ang_dih_imp(nac,ng1,atange,atangm,atdihe,atdihm,
     .                      bondxat,angexat,angmxat,dihexat,dihmxat,
     .               atnamea,nroaa,atxres,atnu,resname,atimp,impxat)

c si existe el residuo WAT, lo pasa a HOH
	do i=1,nac
	 if(aaname(i).eq.'WAT') aaname(i)='HOH'
	enddo

c se fija que las aguas esten en el orden correcto:
	do i=1,nroaa
	 if(resname(i).eq.'HOH') then
	  do j=1,atxres(i)
	   if(atnamea(i,j).eq.'O') then
	    if(j.ne.1) then
	     write(6,*)  'solvent: Wrong order in water residue :',i
	     stop
	    endif
	   endif
	  enddo
	 endif
	enddo

c assign atsinres
        atsinres=0
        if(nroaa.gt.20000) then
        call die('solvent: increase atsinres vector dimension')
        endif
        atsinres(1:nroaa) = atxres(1:nroaa)

        deallocate(atnamea,atnu,resname,atxres)      
        deallocate(qaa,attypea,nataa)
        deallocate(Ema,Rma)

c checking ST and SV parameters
	do i=1,na_u
	if(qmattype(i).ne.'HO'.and.qmattype(i).ne.'HW') then
	if(Rm(i).eq.0.or.Em(i).eq.0) then
	write(6,'(a,i6)') 'solvent: Wrong solute LJ parameter, atom:', i
	STOP
	endif
	endif
	enddo

	do i=1,nac
        if(attype(i).ne.'HO'.and.attype(i).ne.'HW') then
        if(Rm(i+na_u).eq.0.or.Em(i+na_u).eq.0.or.pc(i).eq.0) then
	write(*,*) "Rm(i+na_u)", Rm(i+na_u)
	 write(*,*) "Em(i+na_u)", Em(i+na_u)
	write(*,*) "pc(i)", pc(i)
	write(6,'(a,i6)') 'solvent: Wrong solvent LJ parameter, atom:', i
        STOP
	endif
	endif
	enddo	

       return
 10    stop 'solvent: Problem reading solvent coordinates'
 20    stop 'solvent: Problem reading solute atom types'
 30    stop 'solvent: problem reading cut off radius'
 40    stop 'solvent: Problem reading solvent connectivities'
       end 
c*************************************************************
c  subrutina q pone la carga y tipo de atomo

	subroutine paramats(nroaa,resname,atnamea,qaa,attypea,nataa,
     .                nac,atxres,pc,attype)     
	use ionew
	implicit none
        integer nroaa,nataa(nroaa,100),nac,atxres(nroaa)                                                             
	character*4 atnamea(nroaa,100),resname(nroaa),
     .  attypea(nroaa,100),attype(nac)
	double precision qaa(nroaa,100),pc(0:nac)
        double precision, dimension(:,:,:), allocatable, save ::
     .  pcoord
        double precision, dimension(:,:), allocatable, save ::
     .  pqaa
        character*4 atom
        character*4, dimension(:), allocatable, save::
     .  paanamea
        character*4, dimension(:,:), allocatable, save::
     .  patnamea,pattype
        integer, dimension(:,:), allocatable, save::
     .  atnu,pnataa,patmas
        integer, dimension(:), allocatable, save::
     .  patxres
        integer n1,n2,n3,natoms,pnaas
	integer i,j,k,m
        logical search
        character*12 option 
	integer ui

	integer :: i2,j2,k2 !auxiliars
	logical :: kill
	kill=.false.

	call io_assign(ui) 
        open(unit=ui,file="amber.parm")
 
c lee el nro de atomos, y el nro de residuos
c primer linea de toto tiene nroatomos, nrode residuos
c ahora va a leer todos los  aminoacidos y sus variables
        search=.true.
        do while (search)
        read (ui,*,err=1,end=1) option
        if (option.eq.'residues') then
	read(ui,*,err=1,end=1)  pnaas

        allocate(patnamea(pnaas,100),pcoord(pnaas,100,3),
     .  pattype(pnaas,100),patxres(pnaas),paanamea(pnaas),
     .	patmas(pnaas,100),pqaa(pnaas,100),pnataa(pnaas,100))

	do i=1,pnaas
		read(ui,*,err=1,end=1) paanamea(i), patxres(i)
                do j=1,patxres(i)
                	read(ui,*,err=1,end=1) patnamea(i,j),pattype(i,j),
     .  n1,n2,n3,pnataa(i,j),patmas(i,j),pqaa(i,j)
		enddo
        enddo
        search=.false.
        endif
        enddo
        call io_close(ui)

c si funciona como subrutina asigna las cargas y attypeas seguna amber
	qaa=0.0
        do i=1,nroaa
         do k=1,pnaas
          if(resname(i).eq.paanamea(k)) then
           do j=1,patxres(k)
            do m=1,patxres(k)
             if (atnamea(i,j).eq.patnamea(k,m)) then
              qaa(i,j) = pqaa(k,m)
              attypea(i,j) = pattype(k,m)
              nataa(i,j) = pnataa(k,m)
             endif
            enddo
           enddo
	  else
          endif
         enddo
        enddo

        k=1
        do i=1,nroaa
          do j=1,atxres(i)
            attype(k)=attypea(i,j)
            k=k+1
          enddo
        enddo
 

        k=1
        do i=1,nroaa
          do j=1,atxres(i)
            pc(k)=qaa(i,j)
!		write(6,*) k, attypea(i,j)
	    if(qaa(i,j).eq.0.0) then !report error
	      write(6,'(a,i5)') 'solvent: Wrong atom name  :',k           
	      k2=k

	      if (i.gt.1) then
!		write(*,*) "i>1"
  	        i2=i-1
                do j2=1,atxres(i2)
	          k=k-1
    	        end do
	        i2=i
                do j2=1,j-1
                  k=k-1
                end do


	        do i2=i-1, i
	          do j2=1,atxres(i2)
	            if (i .eq. i2 .and. j2.eq.j) then
		      write(6,*) k, resname(i2), attypea(i2,j2),"<--- this one"
		    else
		      write(6,*) k, resname(i2), attypea(i2,j2)
		    end if
		    k=k+1
		  end do
	        end do

	      else
!		write(*,*) "i=1"
                do j2=1,j-1
                  k=k-1
                end do

                i2=i
                do j2=1,atxres(i2)
                  if (j2.eq.j) then
		    write(6,*) k, resname(i2), attypea(i2,j2),"<--- this one"
                  else
                    write(6,*) k, resname(i2), attypea(i2,j2)
                  end if
                  k=k+1
                end do
	      end if
	      kill=.true.
	      k=k2
	    endif

            k=k+1
          enddo
        enddo


!	write(*,*) "kill?"
	if (kill) STOP

        deallocate(patnamea,pcoord,
     .  pattype,patxres,paanamea,
     .  patmas,pqaa,pnataa) 

        return
 1      stop 
     .'solvent: Problem reading residues block in amber.parm file'
        end
c****************************************************
c subrutina q asigna segun el attypea  los pots Lj 
 
	subroutine lj(nroaa,attypea,Ema,Rma,na_u,natot,
     .                Em,Rm,qmattype,atxres)
	use ionew
	implicit none
        integer i,j,k,ljnum(200),nlj,nroaa,na_u,natot,atxres(nroaa)                                                  
	double precision pEm(200),pRm(200),Ema(nroaa,100),
     .  Rma(nroaa,100),Em(natot),Rm(natot)
        character*4 ljtype(200),attypea(nroaa,100),qmattype(na_u)  
        logical search
        character*12 option 
        integer ui
	Rm=0.d0
	Em=0.d0

        call io_assign(ui)
        open(unit=ui,file="amber.parm")
 
c lee el archivo con los parametros
        search=.true.
        do while (search)
        read (ui,*,err=1,end=1) option
        if (option.eq.'ljs') then
        read(ui,*,err=1,end=1) nlj
        do  i=1,nlj
        if(nlj.ge.200) stop 'solvent: LJ parameters must not exeed 200'
                read (ui,*,err=1,end=1) ljtype(i),pRm(i),pEm(i)
                ljnum(i)=i
        enddo

        search=.false.
        endif
        enddo
        call io_close(ui)

c  pasa los LJ a las unidades del siesta
        do i=1,nlj
        pRm(i) = (2.0d0*pRm(i)/0.529177d0)/(2.d0**(1.d0/6.d0))
        pEm(i) = (pEm(i)/627.5108d0)
        enddo
 
c asigna el LJ corresp al attypea xa el solvente
        do i=1,nroaa
        do j=1,atxres(i)
 
        if (attypea(i,j).eq.'C'.or.
     .  attypea(i,j).eq.'CA'.or.attypea(i,j).eq.'CM'.or.
     .  attypea(i,j).eq.'CC'.or.attypea(i,j).eq.'CV'.or.
     .  attypea(i,j).eq.'CW'.or.attypea(i,j).eq.'CR'.or.
     .  attypea(i,j).eq.'CB'.or.attypea(i,j).eq.'C*'.or.
     .  attypea(i,j).eq.'CN'.or.attypea(i,j).eq.'CK'.or.
     .  attypea(i,j).eq.'CQ'.or.attypea(i,j).eq.'CX'.or.
     .  attypea(i,j).eq.'CY'.or.attypea(i,j).eq.'CD')
     .  then
 
        do k=1,nlj
        if (ljtype(k).eq.'C') then
        Rma(i,j) = pRm(k)
        Ema(i,j) = pEm(k)
        endif
        enddo
 
        elseif (attypea(i,j).eq.'N'.or.attypea(i,j).eq.'NA'.or.
     .  attypea(i,j).eq.'NB'.or.attypea(i,j).eq.'NC'.or.
     .  attypea(i,j).eq.'N*'.or.attypea(i,j).eq.'N2'.or.
     .  attypea(i,j).eq.'NO'.or.attypea(i,j).eq.'NP') then
 
        do k=1,nlj
        if (ljtype(k).eq.'N') then
        Rma(i,j) = pRm(k)
        Ema(i,j) = pEm(k)
        endif
        enddo
 
        else
	do k=1,nlj
	if (attypea(i,j).eq.ljtype(k)) then
 	Rma(i,j) = pRm(k)
        Ema(i,j) = pEm(k)
        endif
        enddo
	
	endif
	enddo
	enddo

c asigna el LJ corresp al qmttype xa el soluto
	do i=1,na_u
 
        if (qmattype(i).eq.'C'.or.
     .  qmattype(i).eq.'CA'.or.qmattype(i).eq.'CM'.or.
     .  qmattype(i).eq.'CC'.or.qmattype(i).eq.'CV'.or.
     .  qmattype(i).eq.'CW'.or.qmattype(i).eq.'CR'.or.
     .  qmattype(i).eq.'CB'.or.qmattype(i).eq.'C*'.or.
     .  qmattype(i).eq.'CN'.or.qmattype(i).eq.'CK'.or.
     .  qmattype(i).eq.'CQ'.or.qmattype(i).eq.'CX'.or.
     .  qmattype(i).eq.'CY'.or.qmattype(i).eq.'CD') then
 
        do k=1,nlj
        if (ljtype(k).eq.'C') then
        Rm(i) = pRm(k)
        Em(i) = pEm(k)
        endif
        enddo

        elseif (qmattype(i).eq.'N' .or.qmattype(i).eq.'NA'.or.
     .          qmattype(i).eq.'NB'.or.qmattype(i).eq.'NC'.or.
     .          qmattype(i).eq.'N*'.or.qmattype(i).eq.'N2'.or.
     .          qmattype(i).eq.'NO'.or.qmattype(i).eq.'NP') then
 
        do k=1,nlj
        if (ljtype(k).eq.'N') then
        Rm(i) = pRm(k)
        Em(i) = pEm(k)
        endif
        enddo
 
        else
        do k=1,nlj
        if (qmattype(i).eq.ljtype(k)) then
        Rm(i) = pRm(k)
        Em(i) = pEm(k)
        endif
        enddo
 
        endif
        enddo

c pasa los LJ del sv
        k=na_u+1
        do i=1,nroaa
        do j=1,atxres(i)
        Em(k)=Ema(i,j)
        Rm(k)=Rma(i,j)     
        k=k+1
        enddo
        enddo

        return
 1      write(*,*)
     .'solvent: Problem reading LJ block in amber.parm file',i
	stop
        end
c*******************************************************************
c subrutina q asigna los 1eros vecinos

       subroutine amber_ng1(nac,nataa,nroaa,atxres,atnu,resname,
     .                      ng1,atnamea,ncon,con)
        use ionew
        use precision
        implicit none
        integer i,j,k,l,m,n,nac,na_u,nresid,nroaa,                                                                   
     .   nataa(nroaa,100),atxres(nroaa),atnu(nroaa,100),                                                             
     .   ng1(nac,6),ncon,con(2,ncon)         
        character*4 resname(nroaa),atnamea(nroaa,100)   
        character*4, dimension(:), allocatable, save::
     .   presname
        integer, dimension(:,:,:), allocatable, save:: 
     .   png1
        integer, dimension(:), allocatable, save::
     .   bondxres
        logical search
        character*12 option
        character c1*1,c4*4,c2*2
        integer ui

        call io_assign(ui)
        open(unit=ui,file="amber.parm")

        search=.true.
        do while (search)
        read (ui,*,err=20,end=20) option
        if (option.eq.'connectivity') then
        read(ui,*,err=20,end=20) nresid
        allocate(presname(nresid),bondxres(nresid),png1(nresid,100,2))

        do i=1,nresid
        read(ui,*,err=20,end=20) presname(i),bondxres(i)
           do j=1,bondxres(i)
           read(ui,*,err=20,end=20) png1(i,j,1),png1(i,j,2)
           enddo
        enddo 

        search=.false.
        endif
        enddo
        call io_close(ui)

c asigna los 1eros vecinos de cada atomo      
         do i=1,nroaa
           do k=1,nresid
             if (resname(i).eq.presname(k)) then
               do j=1,atxres(i)
                 n=1  
                 do l=1,bondxres(k)
                   if (nataa(i,j).eq.png1(k,l,1)) then
                     do m=1,atxres(i)
                       if (nataa(i,m).eq.png1(k,l,2)) then 
                         ng1(atnu(i,j),n) = atnu(i,m)
                         n=n+1
                       endif
                     enddo
                   elseif(nataa(i,j).eq.png1(k,l,2)) then
                     do m=1,atxres(i)
                       if (nataa(i,m).eq.png1(k,l,1)) then
                         ng1(atnu(i,j),n) = atnu(i,m)
                         n=n+1
                       endif
                     enddo 
                   endif 
                 enddo
               enddo
             endif
           enddo
         enddo

c calcula los vecinos entre 2 aa seguidos
c esto es lo q hay q modificar parta omitir bonds entre residuos que no correspondan, Nick
         do i=2,nroaa
         c4=resname(i)
         c1=c4(1:1)
         if(c1.eq.'N') then
         if(c4.eq.'NME') then
         continue
         else 
         goto 5 
         endif
         endif
         do j=1,atxres(i)
         if(atnamea(i,j).eq.'N') then
! se fija si hay exclusion por connectivity extra 1er num -1 segundo el N
                do k=1,ncon
!                write(*,*) k,con(1,k),con(2,k)
                if(con(1,k).eq.-1.and.con(2,k).eq.atnu(i,j)) then
                write(*,*) 'skiping N atom number',atnu(i,j)
                goto 5    
                endif    
                enddo      
         write(*,*) 'Connecting N atom number',atnu(i,j)
         do m=1,atxres(i-1)
         if(atnamea(i-1,m).eq.'C') then
         ng1(atnu(i,j),3) = atnu(i-1,m)
         ng1(atnu(i-1,m),3) = atnu(i,j)
         write(*,*) 'wtih ', atnu(i-1,m)
         endif
         enddo
         endif
         enddo
 5       enddo

c calcula los vecinos entre 2 nucleotidos seguidos
         do i=2,nroaa
         c4=resname(i)
         c1=c4(3:3)
         if(c1.eq.'5') goto 6
         do j=1,atxres(i)
         if(atnamea(i,j).eq.'P') then
         do m=1,atxres(i-1)
         c4=atnamea(i-1,m)
         c2=c4(1:2)
         if(c2.eq.'O3') then
         ng1(atnu(i,j),4) = atnu(i-1,m)
         ng1(atnu(i-1,m),2) = atnu(i,j)
         endif
         enddo
         endif
         enddo
 6       enddo

c asigna las uniones de los atomos impuestos en el imput
        do i=1,ncon
c barre conectivodades extra agregadas 
          if(con(1,i).eq.-1) goto 10
          do k=1,6
            if(ng1(con(1,i),k).eq.0) then
              ng1(con(1,i),k)=con(2,i)
		write(*,*) "agregue conectividad", con(1,i), con(2,i), k
              do j=1,6
                if(ng1(con(2,i),j).eq.0) then        
                  ng1(con(2,i),j)=con(1,i)   
	write(*,*) "agregue conectividad", con(1,i), con(2,i), j
                  goto 10
                endif
              enddo
            endif
          enddo
 10       continue                 
        enddo

        deallocate(presname,bondxres,png1)      
        
        return
 20     stop 
     .'solvent: Problem reading connectivity block in amber.parm file'
        end
c**********************************************************************
c subrutina q calcula los bonds, angles, dihedrals and improper torsions

       subroutine bon_ang_dih_imp(nac,ng1,atange,atangm,atdihe,atdihm,
     .                      bondxat,angexat,angmxat,dihexat,dihmxat,
     .               atnamea,nroaa,atxres,atnu,resname,atimp,impxat)
        use ionew 
        implicit none
	integer nac,nroaa
c       parametros asoc a la asignacion de angulos y dihedros
        integer atange(nac,25,2),atangm(nac,25,2),
     .  atdihe(nac,100,3),atdihm(nac,100,3)
        integer i,j,k,l,m,n,t,t2,bondxat(nac),angexat(nac),
     .  dihexat(nac),dihmxat(nac),angmxat(nac),ng1(nac,6)
c       parametros asoc a la asignacion de impropers
        character*4 resname(nroaa),atnamea(nroaa,100)
        integer  atxres(nroaa),atnu(nroaa,100),impxat(nac),
     .  atimp(nac,25,4),nresid,nataa(nroaa,100),imptot,size
        logical search
        character*10 option
        character*4, dimension(:,:,:), allocatable, save::
     .  impatnamea
        character*4, dimension(:), allocatable, save::
     .  presname
        integer, dimension(:), allocatable, save::
     .  pimpxres
        integer, dimension(:,:), allocatable, save::
     .  impnum
        integer ui

c lee y calcula los impropios
        call io_assign(ui)
        open(unit=ui,file="amber.parm")
        search=.true.
        do while (search)
        read (ui,*,err=1,end=1) option
        if (option.eq.'impropers') then
        read(ui,*,err=1,end=1) nresid
 
        allocate(presname(nresid),pimpxres(nresid),
     .           impatnamea(nresid,50,4))
 
        do i=1,nresid
        read(ui,*,err=1,end=1) presname(i),pimpxres(i)
        do j=1,pimpxres(i)
        read(ui,*,err=1,end=1) impatnamea(i,j,1),impatnamea(i,j,2),
     .  impatnamea(i,j,3),impatnamea(i,j,4)
        enddo
        enddo

        search=.false.
        endif
        enddo
        call io_close(ui)

c asignacion segun el atomo a partir del aa(2)
c asignacion del numero de impropios imxpat
        size=nroaa*25
        imptot=size
        allocate(impnum(imptot,4))
        impnum=0

        imptot=1
        do i=1,nresid
        do m=1,nroaa
        if(presname(i).eq.resname(m)) then
        do j=1,pimpxres(i)
        do k=1,4

        if(impatnamea(i,j,k).eq.'+M'.and.
     .  m.ne.nroaa) then
                do n=1,atxres(m+1)
                if(atnamea(m+1,n).eq.'N') then
                impnum(imptot,k)=atnu(m+1,n)
                endif
                enddo
        elseif(impatnamea(i,j,k).eq.'-M'.and.
     .  m.ne.1) then
                do n=1,atxres(m-1)
                if(atnamea(m-1,n).eq.'C') then
                impnum(imptot,k)=atnu(m-1,n)
                endif
                enddo
        else
                do n=1,atxres(m)
                if(atnamea(m,n).eq.impatnamea(i,j,k)) then
                impnum(imptot,k)=atnu(m,n)
                endif
                enddo
        endif
        enddo
        imptot=imptot+1
        if(imptot.ge.size) then
        stop 'solvent: increase size of improper matrix'
        endif
        enddo
        endif
        enddo
        enddo
 
        do i=1,imptot
        do j=1,4
        if (impnum(i,j).eq.0) then
         do k=1,4
         impnum(i,K)=0
        enddo
        endif
        enddo
        enddo
 
        do i=1,nac
        k=0
        do j=1,imptot
        if (impnum(j,1).eq.i) then
        k=k+1
        atimp(i,k,1)=impnum(j,1)
        atimp(i,k,2)=impnum(j,2)
        atimp(i,k,3)=impnum(j,3)
        atimp(i,k,4)=impnum(j,4)
        elseif (impnum(j,2).eq.i) then
        k=k+1
        atimp(i,k,1)=impnum(j,1)
        atimp(i,k,2)=impnum(j,2)
        atimp(i,k,3)=impnum(j,3)
        atimp(i,k,4)=impnum(j,4)
        elseif (impnum(j,3).eq.i) then
        k=k+1
        atimp(i,k,1)=impnum(j,1)
        atimp(i,k,2)=impnum(j,2)
        atimp(i,k,3)=impnum(j,3)
        atimp(i,k,4)=impnum(j,4)
        elseif (impnum(j,4).eq.i) then
        k=k+1
        atimp(i,k,1)=impnum(j,1)
        atimp(i,k,2)=impnum(j,2)
        atimp(i,k,3)=impnum(j,3)
        atimp(i,k,4)=impnum(j,4)
        endif
        enddo
        impxat(i)=k
        enddo

        deallocate(impnum,presname,pimpxres,impatnamea)

c aca calcula bonds, angulos y dihedros
         do i=1,nac
         bondxat(i)=0
         do j=1,6
         if (ng1(i,j).ne.0) bondxat(i)=bondxat(i)+1
         enddo
         enddo
 
c       busca angulos con i en la esquina(e)
        do i=1,nac
         k=1
         do j=1,bondxat(i)
          t=ng1(i,j)
          do m=1,bondxat(t)
           if(ng1(t,m).ne.i) then
                 atange(i,k,1)=ng1(i,j)
             atange(i,k,2)=ng1(t,m)
             k=k+1
           endif
          enddo
         enddo
        angexat(i)=k-1
      enddo
 
c       busca angulos con i en el medio(m)
        do i=1,nac
         k=1
           do j=1,bondxat(i)
           do m=1,bondxat(i)
            if(ng1(i,m).gt.ng1(i,j)) then
               atangm(i,k,1)=ng1(i,j)
             atangm(i,k,2)=ng1(i,m)
             k=k+1
            endif
         enddo
         enddo
         angmxat(i)=k-1
        enddo

c      busca los atdihedros con i en el extremo(e)
        do i=1,nac
         k=1
         do j=1,angexat(i)
           t=atange(i,j,2)
          do m=1,bondxat(t)
            t2=ng1(t,m)
           if(t2.ne.atange(i,j,1)) then
                 atdihe(i,k,1)=atange(i,j,1)
             atdihe(i,k,2)=t
             atdihe(i,k,3)=t2
             k=k+1
            endif
           enddo
         enddo
        dihexat(i)=k-1
        enddo
 
c       busca los atdihedros con i en el medio(m)
c       primero a partir de angulos con i en la esquina
        do i=1,nac
         k=1
           do j=1,angexat(i)
           t=atange(i,j,1)
                do m=1,bondxat(i)
             t2=ng1(i,m)
                 if(t2.ne.t) then
                atdihm(i,k,1)=t2
             atdihm(i,k,2)=t
             atdihm(i,k,3)=atange(i,j,2)
             k=k+1
           endif
          enddo
         enddo
          dihmxat(i)=k-1
        enddo

        return
 1      write(*,*)
     .'solvent: Problem reading impropers block in amber.parm file',i,j
	stop
        end
c************************************************************************
c subrutina q lee cuantos atomos tiene cada aa

        subroutine atxaa(n,aanamea,atomsxaa)
        use ionew 
        implicit none
        character*4 aanamea(200),nada
        integer atomsxaa(200),i,j,k,l,m,n
        logical search
        character*10 option
        integer ui
        call io_assign(ui)
        open(unit=ui,file="amber.parm")
        search=.true.
        do while (search)
        read (ui,*,err=1,end=1) option
 
        if (option.eq.'residues') then
        read(ui,*,err=1,end=1) n
	if(n.ge.200) then
	stop 'solvent: Number of residues must not exeed 200'
	endif
        do i=1,n
                read(ui,*,err=1,end=1) aanamea(i),atomsxaa(i)
!aanamea es el nombre del residuo
!atomsxaa es la cantidad de atomos q tiene el residuo
!		write(*,*) "lei: ", aanamea(i),atomsxaa(i)
                do j=1,atomsxaa(i)
	if(atomsxaa(i).ge.100) then
	stop 'solvent: Number of atoms in a residue must not exeed 100'
	endif
                read(ui,*,err=1,end=1) nada
                enddo
	enddo
        search=.false.
        endif
        enddo
        call io_close(ui)
        return
 1      write(*,*)
     .'solvent: Problem reading residues block in amber.parm file',i,j
	stop
        end
c**********************************************************************
c subrutina q lee el numero y los parametros de bonds, angles, 
c dihedrals and improper torsions 

        subroutine amber_union_parms(nbond,kbond,bondeq,bondtype,
     .        nangle,kangle,angleeq,angletype,ndihe,kdihe,
     .        diheeq,dihetype,multidihe,perdihe,nimp,kimp,
     .        impeq,imptype,multiimp,perimp,nparm)

        use ionew 
        implicit none
	integer nparm  
        character*5 bondtype(nparm)
        character*8 angletype(nparm)
        character*11 dihetype(nparm),imptype(nparm)
        double precision kbond(nparm),bondeq(nparm),kangle(nparm),
     .  angleeq(nparm),kdihe(nparm),diheeq(nparm),perdihe(nparm),
     .  kimp(nparm),impeq(nparm),perimp(nparm)
        integer nbond,nangle,ndihe,nimp,i,j,k,multidihe(nparm),
     .          multiimp(nparm)
        character option*12
        logical search
        integer ui

c bonds
        call io_assign(ui)
        open(unit=ui,file="amber.parm")
        search=.true.
	bondtype=""
c pongo en 0, Nick
        do while (search)
        read (ui,*,err=100,end=100) option
        if (option.eq.'bonds') then
		read(ui,*,err=100,end=100) nbond
		if(nbond.gt.nparm) stop 'solvent: Increase nparm'
        	do i=1,nbond
        	read(ui,10,err=100,end=100) bondtype(i),kbond(i),bondeq(i)
	        enddo
        search=.false.
        endif
        enddo
        call io_close(ui)

c angles
        call io_assign(ui)
        open(unit=ui,file="amber.parm")
        search=.true.
        do while (search)
        read (ui,*,err=200,end=200) option
        if (option.eq.'angles') then
        	read(ui,*,err=200,end=200) nangle
                if(nangle.gt.nparm) stop 'solvent: Increase nparm'
        	do i=1,nangle
        	read(ui,20,err=200,end=200) angletype(i),kangle(i),angleeq(i)
        	enddo
        search=.false.
        endif
        enddo
        call io_close(ui)

c diherdrals
        call io_assign(ui)
        open(unit=ui,file="amber.parm")
        search=.true.
        do while (search)
        read (ui,*,err=300,end=300) option
        if (option.eq.'dihes') then
        	read(ui,*,err=300,end=300) ndihe
                if(ndihe.gt.nparm) stop 'solvent: Increase nparm'
        	do i=1,ndihe
        	read(ui,30,err=300,end=300) dihetype(i),multidihe(i),kdihe(i),
     .                      diheeq(i),perdihe(i)
        	enddo
        search=.false.
        endif
        enddo
        call io_close(ui)

c impropers
        call io_assign(ui)
        open(unit=ui,file="amber.parm")
        search=.true.
        do while (search)
        read (ui,*,err=400,end=400) option
        if (option.eq.'imps') then
        	read(ui,*,err=400,end=400) nimp
                if(nimp.gt.nparm) stop 'solvent: Increase nparm'
        	do i=1,nimp
        	read(ui,40,err=400,end=400) imptype(i),multiimp(i),kimp(i),
     .                      impeq(i),perimp(i)
	        enddo
        search=.false.
        endif
        enddo
        call io_close(ui)

 10     format(A5,2x,F5.1,4x,F6.4)
 20     format(A8,3x,F5.1,6x,F6.2)
 30     format(A11,3x,I1,3x,F6.3,7x,F7.3,10x,F6.3)
 40     format(A11,3x,I1,3x,F6.3,7x,F7.3,10x,F6.3)
        return
 100    write(*,*) 
     .  'solvent: Problem reading bonds block in amber.parm file',i
        stop
 200    write(*,*) 
     .  'solvent: Problem reading angles block in amber.parm file',i
        stop
 300    write(*,*) 
     .  'solvent: Problem reading dihes block in amber.parm file',i
        stop
 400    write(*,*) 
     .  'solvent: Problem reading imps block in amber.parm file',i
        stop
        end
!**********************************************************************


