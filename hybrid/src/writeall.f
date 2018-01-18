c subroutine that writes the energy  

      subroutine wriene(istp,slabel,idyn,Etots,cfmax)  

      use ionew
      implicit          none
      character         slabel*20, paste*24
      integer           istp,idyn
      double precision  Etots, cfmax
      external          paste 

      character         fname*24 
      logical           frstme
      integer           unit
      double precision  eV,Ang
      save              frstme,Ang,eV,unit,fname
      data              frstme /.true./

c -------------------------------------------------------------------

      if ( frstme ) then
        fname = paste(slabel,'.ene')
        Ang  = 1.d0 / 0.529177d0
        eV  = 1.d0 / 27.211396132d0
c        eV  = 1.d0 / 13.60580d0
c cambiado para E de lio
        frstme = .false.
      endif

      call io_assign(unit)
      open( unit, file=fname, form = 'formatted', position='append',
     .      status='unknown')

	if(idyn .eq. 0 .or. idyn .eq. 5 ) then
	  write(unit,'(i5,2x,F18.7,2x,F14.7)') istp,Etots/eV,cfmax*Ang/eV
	endif

      call io_close(unit)
 
      return
      end

c*************************************************************************
c subroutine that writes the coordinates in PDB/CRD format

      subroutine wripdb(na,slabel,rclas,natot,istp,wricoord,nac,atname,
     .                  aaname,aanum,nesp,atsym,isa,list,block)

      use ionew
      implicit          none
      character         slabel*20, paste*24
      integer           na,nac,natot,istp,wricoord
      integer           nesp, isa(na)
      integer           list(nac),block(nac)
      character*2       atsym(nesp)       
      double precision  rclas(3,natot)
      external          paste

      character         fnamep*24,fnamec*24,fnamel*24,title*30,ch*4
      character*1       chain(natot)
      logical           frstme,foundc
      integer           unitp,unitc,unitl
      double precision  eV,Ang
      save              frstme,Ang,eV,unitc,fnamec,unitl,fnamel
      data              frstme /.true./

c Varibles added to write a pdb file
      integer     i,ia,aanum(nac)
      character*4 aaname(nac),atname(nac)

c -------------------------------------------------------------------
      chain(1:natot)='A'
      do i=1,nac
        if(list(i).eq.0)  chain(i+na)='L'
        if(block(i).eq.0) chain(i+na)='B'
        if(list(i).eq.0.and.block(i).eq.0) chain(i+na)='C'
      enddo

      if ( frstme ) then
        Ang  = 1.d0 / 0.529177d0
        eV  = 1.d0 / 27.211396132d0
!        eV  = 1.d0 / 13.60580d0
      fnamep = paste(slabel,'.init.pdb')
c      fnamec = paste(slabel,'.crd')
      fnamel = paste(slabel,'.last.pdb')

c writes initial PDB coords in .init.pdb file 
      call io_assign(unitp)
      open( unitp, file=fnamep, form = 'formatted', status='unknown')
       write(unitp,'(A4,I7,2x,A2,2x,A4,A,I4,4x,3f8.3)')
     . ('ATOM',i,atsym(isa(i)),'STO ',chain(i),i,
     . (rclas(ia,i)/Ang,ia=1,3), i=1,na)
       write(unitp,'(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)')
     . ('ATOM',na+i,atname(i),aaname(i),chain(na+i),aanum(i),
     . (rclas(ia,na+i)/Ang,ia=1,3), i=1,nac)
       write(unitp,'(A3)') 'END'
       call io_close(unitp)

c writes in .crd file at begining 
c      inquire( file=fnamec, exist=foundc )
c	if (foundc) then
c      call io_assign(unitc)
c      open( unitc, file=fnamec, status='old' )
c      read(unitc,'(a30)',err=1,end=1) title
c      call io_close(unitc)
c      ch=title(1:4)
c      if(ch.ne.'File') then
c      call io_assign(unitc)
c      open( unitc, file=fnamec, form = 'formatted', status='unknown')
c      write(unitc,'(20a)') 'File CRD'
c      call io_close(unitc)
c      endif
c	else
c      call io_assign(unitc)
c      open( unitc, file=fnamec, form = 'formatted', status='unknown')  
c      write(unitc,'(20a)') 'File CRD'
c      call io_close(unitc)
c	endif
c      frstme = .false.
      endif !frstme
c	call flush(6)

c writes in .crd acumulately 
c      if(MOD(istp,wricoord).eq.0) then         
c      call io_assign(unitc)
c      open( unitc, file=fnamec, form = 'formatted', position='append',
c     .      status='unknown')
c       write(unitc,'(10F8.3)') (rclas(1:3,i)/Ang,i=1,natot)
c      call io_close(unitc)
c      endif

c writes actual coords in PDB format in .last.pdb file 
      call io_assign(unitl)
      open( unitl, file=fnamel, form = 'formatted', status='unknown')
       write(unitl,'(A4,I7,2x,A2,2x,A4,A,I4,4x,3f8.3)')
     . ('ATOM',i,atsym(isa(i)),'STO ',chain(i),i,
     . (rclas(ia,i)/Ang,ia=1,3), i=1,na)
       write(unitl,'(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)')
     . ('ATOM',na+i,atname(i),aaname(i),chain(na+i),aanum(i),
     . (rclas(ia,na+i)/Ang,ia=1,3), i=1,nac)
       write(unitl,'(A3,i5)') 'END',istp
       call io_close(unitl)

	return
 1      stop 'write: problem reading from CRD file'
	end

c******************************************************************************
c subrutine that writes the reaction coordinate and its PDB

       subroutine wrirtc(slabel,Etots,rtot,constrpaso,na,nac,natot,
     .                   rclas,atname,aaname,aanum,nesp,atsym,isa)

      use ionew
      implicit          none
      character         slabel*20, paste*24
      integer           na,nac,natot,istp,wricoord
      integer           nesp,isa(na)
      double precision  rclas(3,natot)
      character*2       atsym(nesp)
      double precision  Etots, rtot   
      external          paste
 
      character         fnamee*24,fnamec*24 
      logical           frstme
      integer           i,j,ia,unite,unitc
      double precision  eV, Ang
      save              frstme,Ang,eV,unite,fnamee,unitc,fnamec
      data              frstme /.true./
      integer           aanum(nac),constrpaso
      character*4       aaname(nac),atname(nac)

c -------------------------------------------------------------------

      if ( frstme ) then
        Ang  = 1.d0 / 0.529177d0
        eV  = 1.d0 / 27.211396132d0
      fnamec = paste(slabel,'.rcg')
      fnamee = paste(slabel,'.rce')
      frstme = .false.
      endif

c writes .rce file
      call io_assign(unite)
      open( unite, file=fnamee, form = 'formatted', position='append',
     .      status='unknown')
       write(unite,'(F10.4,2x,F14.7)') rtot, Etots/eV
      call io_close(unite)

c wirtes .rcg file
      call io_assign(unitc)
      open( unitc, file=fnamec, form = 'formatted', position='append',
     .      status='unknown')
       write(unitc,'(A4,I7,2x,A2,2x,A4,A,I4,4x,3f8.3)')
     . ('ATOM',i,atsym(isa(i)),'STO ','A',i,
     . (rclas(ia,i)/Ang,ia=1,3), i=1,na)
       write(unitc,'(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)')
     . ('ATOM',na+i,atname(i),aaname(i),'A',aanum(i),
     . (rclas(ia,na+i)/Ang,ia=1,3), i=1,nac)
      write(unitc,'(A3,i3)') 'END',constrpaso
      call io_close(unitc) 

      return
      end

c***********************************************************************************
	subroutine wrtcrd(natot,rclas)

	use ionew
	use fdf    
	use sys
	implicit none
	integer natot
	integer i,unit,iunit,nconstr,iconstr,typeconstr(20)
	integer atmsconstr(20,20),ndists(20)
	double precision rclas(3,natot),coef(20,10),rtot(20)
	character exp
        character slabel*24, paste*24,fname*24
        external  paste
        logical   frstme,constrlog
        data      frstme /.true./
        data      constrlog /.true./
        integer at1,at2,at3,at4,at5,at6,at7,at8
        double precision rp(3),r12,r34,r56,r78,dist,angle,dihedro2
        save frstme,fname,atmsconstr,ndists,coef,nconstr,typeconstr,
     .       constrlog

c 8 ways of defining the reaction coordinate 
c 1 = r1 - r2 coupled
c 2 = distance 
c 3 = angle
c 4 = dihedral
c 5 = r1 + r2 coupled
c 6 = ( r1 + r2 ) - ( r3 + r4 ) coupled
c 7 = plane to atom distance
c 8 = c1*r1 + c2*r2 + c3*r3 + ....

	if(frstme) then
c read variables
	if ( fdf_block('WrtCrd',iunit) ) then
	read(iunit,'(A)',advance='no',err=100,end=100) exp
	if(exp.eq.'%') then 
	constrlog=.false.
	goto 10
	endif
	read(iunit,*,err=100,end=100) exp,nconstr
	if(nconstr.gt.20) then
	call die('wrtcrd: nconstr must be lower than 20')
	endif
	do iconstr=1,nconstr
	read(iunit,*,err=100,end=100) exp,typeconstr(iconstr)

         if     (typeconstr(iconstr).eq.1) then          
         read(iunit,*,err=100,end=100) exp,(atmsconstr(iconstr,i),i=1,4)
         elseif (typeconstr(iconstr).eq.2) then
         read(iunit,*,err=100,end=100) exp,(atmsconstr(iconstr,i),i=1,2)
         elseif (typeconstr(iconstr).eq.3) then
         read(iunit,*,err=100,end=100) exp,(atmsconstr(iconstr,i),i=1,3)
         elseif (typeconstr(iconstr).eq.4) then
         read(iunit,*,err=100,end=100) exp,(atmsconstr(iconstr,i),i=1,4)
         elseif (typeconstr(iconstr).eq.5) then
         read(iunit,*,err=100,end=100) exp,(atmsconstr(iconstr,i),i=1,4)
         elseif (typeconstr(iconstr).eq.6) then
         read(iunit,*,err=100,end=100) exp,(atmsconstr(iconstr,i),i=1,8) 
         elseif (typeconstr(iconstr).eq.7) then
         read(iunit,*,err=100,end=100) exp,(atmsconstr(iconstr,i),i=1,5)
         elseif (typeconstr(iconstr).eq.8) then
         read(iunit,*,err=100,end=100) exp,ndists(iconstr)

       if(ndists(iconstr).gt.10) then
       call die('wrtcrd: ndists in typeconstr 8 must not exceed 10')
       endif
       read(iunit,*,err=100,end=100) 
     . exp,(coef(iconstr,i),i=1,ndists(iconstr))
       read(iunit,*,err=100,end=100) 
     . exp,(atmsconstr(iconstr,i),i=1,ndists(iconstr)*2)
         else
         call die('wrtcrd: typeconstr must be 1-8')
         endif

	if(i.gt.20) then
	call die('wrtcrd: atoms with constrain must be lower than 20')
	endif

	enddo !nconstr
	else
	constrlog=.false.
	goto 10
	endif !fdf

c name file
        slabel = fdf_string( 'SystemLabel', 'siesta' )
        fname = paste(slabel,'.wrt')

10	continue
	frstme=.false.
	endif !frstme

c write only if constrlog is true
        if(constrlog) then

c change units
        rclas(1:3,1:natot)=rclas(1:3,1:natot)*0.529177d0

c loop over nconstr
        do iconstr=1,nconstr

        if (typeconstr(iconstr).eq.1) then
        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
        at3=atmsconstr(iconstr,3)
        at4=atmsconstr(iconstr,4)  
 
        r12=dist(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .   rclas(1,at2),rclas(2,at2),rclas(3,at2))
        r34=dist(rclas(1,at3),rclas(2,at3),rclas(3,at3),
     .   rclas(1,at4),rclas(2,at4),rclas(3,at4))
        rtot(iconstr)=r34-r12

        elseif (typeconstr(iconstr).eq.2) then
        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
 
        rtot(iconstr)=dist(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .   rclas(1,at2),rclas(2,at2),rclas(3,at2))

        elseif(typeconstr(iconstr).eq.3) then
        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
        at3=atmsconstr(iconstr,3)

       rtot(iconstr)=angle(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .            rclas(1,at2),rclas(2,at2),rclas(3,at2),
     .            rclas(1,at3),rclas(2,at3),rclas(3,at3))

        elseif(typeconstr(iconstr).eq.4) then
        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
        at3=atmsconstr(iconstr,3)
        at4=atmsconstr(iconstr,4)

        rtot(iconstr)=dihedro2(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .   rclas(1,at2),rclas(2,at2),
     .   rclas(3,at2),
     .   rclas(1,at3),rclas(2,at3),
     .   rclas(3,at3),
     .   rclas(1,at4),rclas(2,at4),
     .   rclas(3,at4))

        elseif (typeconstr(iconstr).eq.5) then
        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
        at3=atmsconstr(iconstr,3)
        at4=atmsconstr(iconstr,4)

        r12=dist(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .   rclas(1,at2),rclas(2,at2),rclas(3,at2))
        r34=dist(rclas(1,at3),rclas(2,at3),rclas(3,at3),
     .   rclas(1,at4),rclas(2,at4),rclas(3,at4))
        rtot(iconstr)=r34+r12

        elseif (typeconstr(iconstr).eq.6) then
        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
        at3=atmsconstr(iconstr,3)
        at4=atmsconstr(iconstr,4)
        at5=atmsconstr(iconstr,5)
        at6=atmsconstr(iconstr,6)
        at7=atmsconstr(iconstr,7)
        at8=atmsconstr(iconstr,8)

        r12=dist(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .   rclas(1,at2),rclas(2,at2),rclas(3,at2))
        r34=dist(rclas(1,at3),rclas(2,at3),rclas(3,at3),
     .   rclas(1,at4),rclas(2,at4),rclas(3,at4))
        r56=dist(rclas(1,at5),rclas(2,at5),rclas(3,at5),
     .   rclas(1,at6),rclas(2,at6),rclas(3,at6))
        r78=dist(rclas(1,at7),rclas(2,at7),rclas(3,at7),
     .   rclas(1,at8),rclas(2,at8),rclas(3,at8))
        rtot(iconstr)=(r34+r12)-(r56+r78)

        elseif (typeconstr(iconstr).eq.7) then
        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
        at3=atmsconstr(iconstr,3)
        at4=atmsconstr(iconstr,4)
        at5=atmsconstr(iconstr,5)

        rp(1)=(rclas(1,at2)+rclas(1,at3)+rclas(1,at4)+rclas(1,at5))/4.0
        rp(2)=(rclas(2,at2)+rclas(2,at3)+rclas(2,at4)+rclas(2,at5))/4.0
        rp(3)=(rclas(3,at2)+rclas(3,at3)+rclas(3,at4)+rclas(3,at5))/4.0

        rtot(iconstr)=dist(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .   rp(1),rp(2),rp(3))

        elseif (typeconstr(iconstr).eq.8) then
        rtot(iconstr)=0.0
        do i=1,ndists(iconstr)
        at1=atmsconstr(iconstr,i)
        at2=atmsconstr(iconstr,i+1)

        rtot(iconstr)=dist(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .  rclas(1,at2),rclas(2,at2),rclas(3,at2))

        rtot(iconstr)=rtot(iconstr)+coef(iconstr,i)*rtot(iconstr)

        enddo !ndists
        endif !typeconstr
        enddo  !nconstr

c write file
        call io_assign(unit)
        open( unit, file=fname, form = 'formatted', position='append',
     .        status='unknown')
        write(unit,'(20f14.8)') rtot(1:nconstr)
        call io_close(unit)

c change units 
        rclas(1:3,1:natot)=rclas(1:3,1:natot)/0.529177d0

        endif !constrlog
        return
 100    stop 'wrtcrd: problem reading WrtCrd block'
        end

c************************************************************************************************


        subroutine write_xyz(natot, na_u, iza, pc, rclas)
        use precision, only: dp
        implicit none
        integer, intent(in) :: natot, na_u
        integer, intent(in), dimension(na_u) :: iza
        double precision, intent(in), dimension(natot-na_u) :: pc
        double precision, intent(in), dimension(3,natot) :: rclas
        integer :: inick
        real(dp) :: Ang
        Ang    = 1._dp / 0.529177_dp
        write(34,*) natot
        write(34,*)
        do inick=1,natot
          if (inick.le.na_u) then
            write(34,345) iza(inick), rclas(1:3,inick)/Ang
          else
            write(34,346) pc(inick-na_u), rclas(1:3,inick)/Ang
          end if
        enddo
 345  format(2x, I2,    2x, 3(f10.6,2x))
 346  format(2x, f10.6, 2x, 3(f10.6,2x))
        end subroutine write_xyz
