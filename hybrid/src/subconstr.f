	subroutine subconstr1(nconstr,typeconstr,kforce,nstepconstr,
     .             rini,rfin,atmsconstr,dr,ro,ndists,coef,constrlog)

	use ionew
	use fdf    
	use sys
	implicit none
	integer i,unit,iunit,nconstr,iconstr,typeconstr(20),k
	integer nstepconstr,atmsconstr(20,20),ndists(20)
	double precision kforce(20),ro(20),rini,rfin,dr,zo,coef(20,10)
	character exp
        character slabel*24, paste*24,fname*24
        external  paste
        logical   found, constrlog
        data      found /.false./

c 8 ways of defining the reaction coordinate 
c 1 = r2 - r1 coupled
c 2 = distance 
c 3 = angle
c 4 = dihedral
c 5 = r1 + r2 coupled
c 6 = ( r1 + r2 ) - ( r3 + r4 ) coupled
c 7 = plane to atom distance
c 8 = c1*r1 + c2*r2 + c3*r3 + ....

c read variables
	if ( fdf_block('ConstrainedOpt',iunit) ) then
	read(iunit,'(A)',advance='no',err=100,end=100) exp
	if(exp.eq.'%') then
	constrlog=.false.
	nstepconstr=0
	return
	endif
	read(iunit,*,err=100,end=100) exp,nconstr
	if(nconstr.gt.20) then
	call die('constr opt: nconstr must be lower than 20')
	endif
	read(iunit,*,err=100,end=100) exp,nstepconstr
	if(nstepconstr.eq.0) then
	call die('constr opt: nstepconstr must be larger than 0')
	elseif(nstepconstr.gt.100) then
	call die('constr opt: nstepconstr must be lower than 100')
	endif
	do iconstr=1,nconstr
	read(iunit,*,err=100,end=100) exp,typeconstr(iconstr)
	read(iunit,*,err=100,end=100) exp,kforce(iconstr)
	if(iconstr.eq.1) then
	read(iunit,*,err=100,end=100) exp,rini,exp,rfin
	else
	read(iunit,*,err=100,end=100) exp,ro(iconstr)
	endif

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
       call die('constr opt: ndists in typeconstr 8 must not exceed 10')
       endif
       read(iunit,*,err=100,end=100) 
     . exp,(coef(iconstr,i),i=1,ndists(iconstr))
       read(iunit,*,err=100,end=100) 
     . exp,(atmsconstr(iconstr,i),i=1,ndists(iconstr)*2)
         else
         call die('constr opt: typeconstr must be 1-8')
         endif

	if(i.gt.20) then
	call die('constr opt: atoms with constrain must be lower than 20')
	endif

	enddo
	else
	call die('constr opt: You must specify the ConstraindeOpt block')
	endif
	write(6,'(/,a)') 'constr opt: Starting a Constrained 
     . optimization run'

c calculates initial ro for all types of constraints
        dr=(rfin-rini)/nstepconstr
        ro(1)=rini
        if(rfin.eq.rini) nstepconstr=0

c reads from .rce of a former run
        slabel = fdf_string( 'SystemLabel', 'siesta' )
        fname = paste(slabel,'.rce')
        inquire( file=fname, exist=found )
        if(found) then
        call io_assign(unit)
        open( unit, file=fname )
        k=0
 10     continue
        k=k+1
        read(unit,*,err=200,end=20) zo
        goto 10
 20     continue
        call io_close(unit)
        k=k-1
        if(nstepconstr.gt.k) then
        nstepconstr=nstepconstr-k
        ro(1)=ro(1)+dr*k  
        write(6,'(/,a)') 
     .  'constr opt: Re-starting a Constrained optimization run'
        endif
        endif !found

        return
 100    stop 'constr opt: problem reading ConstrainedOpt block'
 200    stop 'constr opt: problem reading form rce file'
        end

c****************************************************************************
	subroutine subconstr2(nconstr,typeconstr,kforce,rini,rfin,ro,rt,
     .  nstepconstr,atmsconstr,natot,rclas,fdummy,istep,istepconstr,
     .  ndists,coef)

	use sys
	implicit none
	integer i,natot,iconstr,istep,istepconstr,ndists(20)
	integer nconstr,typeconstr(20),atmsconstr(20,20),nstepconstr
	double precision kforce(20),ro(20),rt(20),rini,rfin,coef(20,10)
	
	integer at1,at2,at3,at4,at5,at6,at7,at8
	double precision rclas(3,natot),fdummy(3,natot)
	double precision fce,fnew(3,10),rp(3),r12,r34,r56,r78,dist,
     .  dx,dy,dz
	double precision pi,fdihe(12),F,dihedro2,rtot,req,kf
	double precision angle,scal,scalar,dscalar,r32,dr12r32
	pi = acos(-1.0d0)

c change units 
        rclas(1:3,1:natot)=rclas(1:3,1:natot)*0.529177d0

c loop over nconstr
        do iconstr=1,nconstr
        fnew = 0.d0

c loop over typeconstr
      if (typeconstr(iconstr).eq.1) then

        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
        at3=atmsconstr(iconstr,3)
        at4=atmsconstr(iconstr,4)  
 
        r12=dist(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .   rclas(1,at2),rclas(2,at2),rclas(3,at2))
        r34=dist(rclas(1,at3),rclas(2,at3),rclas(3,at3),
     .   rclas(1,at4),rclas(2,at4),rclas(3,at4))
        rtot=r34-r12

        rt(iconstr)=rtot
        req=ro(iconstr) 
        kf=kforce(iconstr)

c atom1: dr12
        dx=0.d0
        dy=0.d0
        dz=0.d0
          dx=(1.d0/r12)*(rclas(1,at1)-rclas(1,at2))
          dx=2.d0*kf*(rtot-req)*dx
          dy=(1.d0/r12)*(rclas(2,at1)-rclas(2,at2))
          dy=2.d0*kf*(rtot-req)*dy
          dz=(1.d0/r12)*(rclas(3,at1)-rclas(3,at2))
          dz=2.d0*kf*(rtot-req)*dz
        fnew(1,1)=dx
        fnew(2,1)=dy
        fnew(3,1)=dz

c atom3: dr34
        dx=0.d0
        dy=0.d0
        dz=0.d0
          dx=(1.d0/r34)*(rclas(1,at3)-rclas(1,at4))
          dx=2.d0*kf*(rtot-req)*dx
          dy=(1.d0/r34)*(rclas(2,at3)-rclas(2,at4))
          dy=2.d0*kf*(rtot-req)*dy
          dz=(1.d0/r34)*(rclas(3,at3)-rclas(3,at4))
          dz=2.d0*kf*(rtot-req)*dz
        fnew(1,3)=-dx
        fnew(2,3)=-dy
        fnew(3,3)=-dz

c atom2: -dr12
        fnew(1,2)=-fnew(1,1)
        fnew(2,2)=-fnew(2,1)
        fnew(3,2)=-fnew(3,1)

c atom4: -dr34
        fnew(1,4)=-fnew(1,3)
        fnew(2,4)=-fnew(2,3)
        fnew(3,4)=-fnew(3,3)

c adding fnew to fdummy 
        fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)
        fdummy(1:3,at2)= fdummy(1:3,at2)+fnew(1:3,2)
        fdummy(1:3,at3)= fdummy(1:3,at3)+fnew(1:3,3) 
        fdummy(1:3,at4)= fdummy(1:3,at4)+fnew(1:3,4)      

c test write de fuerzas
	write(666,*) "r, Frestr", rtot, -kf*(rtot-req)



      elseif (typeconstr(iconstr).eq.2) then

        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
 
        rtot=dist(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .   rclas(1,at2),rclas(2,at2),rclas(3,at2))

        rt(iconstr)=rtot
        req=ro(iconstr)
        kf=kforce(iconstr)
 
c atom1: dr12
        dx=0.d0
        dy=0.d0
        dz=0.d0
          dx=(1.d0/rtot)*(rclas(1,at1)-rclas(1,at2))
          dx=2.d0*kf*(rtot-req)*dx
          dy=(1.d0/rtot)*(rclas(2,at1)-rclas(2,at2))
          dy=2.d0*kf*(rtot-req)*dy
          dz=(1.d0/rtot)*(rclas(3,at1)-rclas(3,at2))
          dz=2.d0*kf*(rtot-req)*dz
        fnew(1,1)=-dx
        fnew(2,1)=-dy
        fnew(3,1)=-dz

c atom2: -fce over atom1
        fnew(1,2)=-fnew(1,1)
        fnew(2,2)=-fnew(2,1)
        fnew(3,2)=-fnew(3,1)
 
c adding fnew to fdummy 
        fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)
        fdummy(1:3,at2)= fdummy(1:3,at2)+fnew(1:3,2)


        write(666,*) "r, Frestr", rtot, -kf*(rtot-req)

      elseif(typeconstr(iconstr).eq.3) then
        
        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
        at3=atmsconstr(iconstr,3)

        rtot =angle(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .            rclas(1,at2),rclas(2,at2),rclas(3,at2),
     .            rclas(1,at3),rclas(2,at3),rclas(3,at3))

        rt(iconstr)=rtot
        req=ro(iconstr)
        kf=kforce(iconstr)

       scal=scalar(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .             rclas(1,at2),rclas(2,at2),rclas(3,at2),
     .             rclas(1,at3),rclas(2,at3),rclas(3,at3))
       r12=dist(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .          rclas(1,at2),rclas(2,at2),rclas(3,at2))
       r32=dist(rclas(1,at3),rclas(2,at3),rclas(3,at3),
     .          rclas(1,at2),rclas(2,at2),rclas(3,at2))
       fce=2.d0*kf*(rtot-req)*pi/180d0

c atom1:
       dscalar=(rclas(1,at3)-rclas(1,at2))
       dr12r32=r32*(rclas(1,at1)-rclas(1,at2))/(r12)
       dx=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
       dx=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dx
       dx=fce*dx
       dscalar=(rclas(2,at3)-rclas(2,at2))
       dr12r32=r32*(rclas(2,at1)-rclas(2,at2))/(r12)
       dy=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
       dy=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dy
       dy=fce*dy
       dscalar=(rclas(3,at3)-rclas(3,at2))
       dr12r32=r32*(rclas(3,at1)-rclas(3,at2))/(r12)
       dz=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
       dz=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dz
       dz=fce*dz
       fnew(1,1)=-dx
       fnew(2,1)=-dy
       fnew(3,1)=-dz

c atom2:
       dscalar=2.d0*rclas(1,at2)-rclas(1,at1)-rclas(1,at3)
       dr12r32=(r32*(-rclas(1,at1)+rclas(1,at2))/r12)+
     .         (r12*(-rclas(1,at3)+rclas(1,at2))/r32)
       dx=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
       dx=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dx
       dx=fce*dx
       dscalar=2.d0*rclas(2,at2)-rclas(2,at1)-rclas(2,at3)
       dr12r32=(r32*(-rclas(2,at1)+rclas(2,at2))/r12)+
     .         (r12*(-rclas(2,at3)+rclas(2,at2))/r32)
       dy=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
       dy=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dy
       dy=fce*dy
       dscalar=2.d0*rclas(3,at2)-rclas(3,at1)-rclas(3,at3)
       dr12r32=(r32*(-rclas(3,at1)+rclas(3,at2))/r12)+
     .         (r12*(-rclas(3,at3)+rclas(3,at2))/r32)
       dz=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
       dz=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dz
       dz=fce*dz
       fnew(1,2)=-dx
       fnew(2,2)=-dy
       fnew(3,2)=-dz

c atom3:
       fnew(1,3)=-fnew(1,1)
       fnew(2,3)=-fnew(2,1)
       fnew(3,3)=-fnew(3,1)

c adding fnew to fdummy
        fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)
        fdummy(1:3,at2)= fdummy(1:3,at2)+fnew(1:3,2)
        fdummy(1:3,at3)= fdummy(1:3,at3)+fnew(1:3,3)

      elseif(typeconstr(iconstr).eq.4) then

        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
        at3=atmsconstr(iconstr,3)
        at4=atmsconstr(iconstr,4)

        rtot=dihedro2(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .   rclas(1,at2),rclas(2,at2),
     .   rclas(3,at2),
     .   rclas(1,at3),rclas(2,at3),
     .   rclas(3,at3),
     .   rclas(1,at4),rclas(2,at4),
     .   rclas(3,at4))

        req=ro(iconstr)
        kf=kforce(iconstr)

c forces
        if(req.lt.90 .and.rtot.gt.180) rtot=rtot-360.d0
        if(req.gt.270.and.rtot.lt.180) rtot=rtot+360.d0 
        F=2.d0*kf*(rtot-req)*pi/180.d0
        call diheforce2(natot,rclas,at1,at2,at3,at4,
     .          1,F,fdihe)
        fnew(1:3,1)=fdihe(1:3)
        call diheforce2(natot,rclas,at1,at2,at3,at4,
     .          2,F,fdihe)
        fnew(1:3,2)=fdihe(4:6)
        call diheforce2(natot,rclas,at1,at2,at3,at4,
     .          3,F,fdihe)
        fnew(1:3,3)=fdihe(7:9)
        call diheforce2(natot,rclas,at1,at2,at3,at4,
     .          4,F,fdihe)
        fnew(1:3,4)=fdihe(10:12)

c adding fnew to fdummy
        if((rtot.ge.0..and.rtot.le.180.).or.(rtot.gt.360)) then
         fnew(1:3,1:4)=(-1.d0)*fnew(1:3,1:4)
        elseif((rtot.gt.180..and.rtot.lt.360).or.(rtot.lt.0)) then
         fnew(1:3,1:4)=fnew(1:3,1:4)
        else
         call die('constr opt: Wrong dihedral angle value')
        endif
        fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)
        fdummy(1:3,at2)= fdummy(1:3,at2)+fnew(1:3,2)
        fdummy(1:3,at3)= fdummy(1:3,at3)+fnew(1:3,3)
        fdummy(1:3,at4)= fdummy(1:3,at4)+fnew(1:3,4)
        rt(iconstr)=rtot

      elseif (typeconstr(iconstr).eq.5) then

        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
        at3=atmsconstr(iconstr,3)
        at4=atmsconstr(iconstr,4)

        r12=dist(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .   rclas(1,at2),rclas(2,at2),rclas(3,at2))
        r34=dist(rclas(1,at3),rclas(2,at3),rclas(3,at3),
     .   rclas(1,at4),rclas(2,at4),rclas(3,at4))
        rtot=r34+r12

        rt(iconstr)=rtot
        req=ro(iconstr)
        kf=kforce(iconstr)

c atom1: dr12
        dx=0.d0
        dy=0.d0
        dz=0.d0
          dx=(1.d0/r12)*(rclas(1,at1)-rclas(1,at2))
          dx=2.d0*kf*(rtot-req)*dx
          dy=(1.d0/r12)*(rclas(2,at1)-rclas(2,at2))
          dy=2.d0*kf*(rtot-req)*dy
          dz=(1.d0/r12)*(rclas(3,at1)-rclas(3,at2))
          dz=2.d0*kf*(rtot-req)*dz
        fnew(1,1)=-dx
        fnew(2,1)=-dy
        fnew(3,1)=-dz

c atom3: dr34
        dx=0.d0
        dy=0.d0
        dz=0.d0
          dx=(1.d0/r34)*(rclas(1,at3)-rclas(1,at4))
          dx=2.d0*kf*(rtot-req)*dx
          dy=(1.d0/r34)*(rclas(2,at3)-rclas(2,at4))
          dy=2.d0*kf*(rtot-req)*dy
          dz=(1.d0/r34)*(rclas(3,at3)-rclas(3,at4))
          dz=2.d0*kf*(rtot-req)*dz
        fnew(1,3)=-dx
        fnew(2,3)=-dy
        fnew(3,3)=-dz

c atom2: -dr12
        fnew(1,2)=-fnew(1,1)
        fnew(2,2)=-fnew(2,1)
        fnew(3,2)=-fnew(3,1)

c atom4: -dr34
        fnew(1,4)=-fnew(1,3)
        fnew(2,4)=-fnew(2,3)
        fnew(3,4)=-fnew(3,3)

c adding fnew to fdummy
        fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)
        fdummy(1:3,at2)= fdummy(1:3,at2)+fnew(1:3,2)
        fdummy(1:3,at3)= fdummy(1:3,at3)+fnew(1:3,3)
        fdummy(1:3,at4)= fdummy(1:3,at4)+fnew(1:3,4)

c test write de fuerzas
        write(666,*) "r, Frestr", rtot, -kf*(rtot-req)

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
        rtot=(r34+r12)-(r56+r78)

        rt(iconstr)=rtot
        req=ro(iconstr)
        kf=kforce(iconstr)

c atom1: dr12
        dx=0.d0
        dy=0.d0
        dz=0.d0
          dx=(1.d0/r12)*(rclas(1,at1)-rclas(1,at2))
          dx=2.d0*kf*(rtot-req)*dx
          dy=(1.d0/r12)*(rclas(2,at1)-rclas(2,at2))
          dy=2.d0*kf*(rtot-req)*dy
          dz=(1.d0/r12)*(rclas(3,at1)-rclas(3,at2))
          dz=2.d0*kf*(rtot-req)*dz
        fnew(1,1)=-dx
        fnew(2,1)=-dy
        fnew(3,1)=-dz

c atom3: dr34
        dx=0.d0
        dy=0.d0
        dz=0.d0
          dx=(1.d0/r34)*(rclas(1,at3)-rclas(1,at4))
          dx=2.d0*kf*(rtot-req)*dx
          dy=(1.d0/r34)*(rclas(2,at3)-rclas(2,at4))
          dy=2.d0*kf*(rtot-req)*dy
          dz=(1.d0/r34)*(rclas(3,at3)-rclas(3,at4))
          dz=2.d0*kf*(rtot-req)*dz
        fnew(1,3)=-dx
        fnew(2,3)=-dy
        fnew(3,3)=-dz

c atom2: -dr12
        fnew(1,2)=-fnew(1,1)
        fnew(2,2)=-fnew(2,1)
        fnew(3,2)=-fnew(3,1)

c atom4: -dr34
        fnew(1,4)=-fnew(1,3)
        fnew(2,4)=-fnew(2,3)
        fnew(3,4)=-fnew(3,3)

c atom5: dr56
        dx=0.d0
        dy=0.d0
        dz=0.d0
          dx=(1.d0/r56)*(rclas(1,at5)-rclas(1,at6))
          dx=2.d0*kf*(rtot-req)*dx
          dy=(1.d0/r56)*(rclas(2,at5)-rclas(2,at6))
          dy=2.d0*kf*(rtot-req)*dy
          dz=(1.d0/r56)*(rclas(3,at5)-rclas(3,at6))
          dz=2.d0*kf*(rtot-req)*dz
        fnew(1,5)=dx
        fnew(2,5)=dy
        fnew(3,5)=dz

c atom7: dr78
        dx=0.d0
        dy=0.d0
        dz=0.d0
          dx=(1.d0/r78)*(rclas(1,at7)-rclas(1,at8))
          dx=2.d0*kf*(rtot-req)*dx
          dy=(1.d0/r78)*(rclas(2,at7)-rclas(2,at8))
          dy=2.d0*kf*(rtot-req)*dy
          dz=(1.d0/r78)*(rclas(3,at7)-rclas(3,at8))
          dz=2.d0*kf*(rtot-req)*dz
        fnew(1,7)=dx
        fnew(2,7)=dy
        fnew(3,7)=dz

c atom6: -dr56
        fnew(1,6)=-fnew(1,5)
        fnew(2,6)=-fnew(2,5)
        fnew(3,6)=-fnew(3,5)

c atom8: -dr78
        fnew(1,8)=-fnew(1,7)
        fnew(2,8)=-fnew(2,7)
        fnew(3,8)=-fnew(3,7)

c adding fnew to fdummy
        fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)
        fdummy(1:3,at2)= fdummy(1:3,at2)+fnew(1:3,2)
        fdummy(1:3,at3)= fdummy(1:3,at3)+fnew(1:3,3)
        fdummy(1:3,at4)= fdummy(1:3,at4)+fnew(1:3,4)
        fdummy(1:3,at5)= fdummy(1:3,at5)+fnew(1:3,5)
        fdummy(1:3,at6)= fdummy(1:3,at6)+fnew(1:3,6)
        fdummy(1:3,at7)= fdummy(1:3,at7)+fnew(1:3,7)
        fdummy(1:3,at8)= fdummy(1:3,at8)+fnew(1:3,8)

      elseif (typeconstr(iconstr).eq.7) then

        at1=atmsconstr(iconstr,1)
        at2=atmsconstr(iconstr,2)
        at3=atmsconstr(iconstr,3)
        at4=atmsconstr(iconstr,4)
        at5=atmsconstr(iconstr,5)

        rp(1)=(rclas(1,at2)+rclas(1,at3)+rclas(1,at4)+rclas(1,at5))/4.d0
        rp(2)=(rclas(2,at2)+rclas(2,at3)+rclas(2,at4)+rclas(2,at5))/4.d0
        rp(3)=(rclas(3,at2)+rclas(3,at3)+rclas(3,at4)+rclas(3,at5))/4.d0

        rtot=dist(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .   rp(1),rp(2),rp(3))

        rt(iconstr)=rtot
        req=ro(iconstr)
        kf=kforce(iconstr)

c atom1: dr12
        dx=0.d0
        dy=0.d0
        dz=0.d0
          dx=(1.d0/rtot)*(rclas(1,at1)-rp(1))
          dx=2.d0*kf*(rtot-req)*dx
          dy=(1.d0/rtot)*(rclas(2,at1)-rp(2))
          dy=2.d0*kf*(rtot-req)*dy
          dz=(1.d0/rtot)*(rclas(3,at1)-rp(3))
          dz=2.d0*kf*(rtot-req)*dz
        fnew(1,1)=-dx
        fnew(2,1)=-dy
        fnew(3,1)=-dz

c adding fnew to fdummy
        fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)

      elseif (typeconstr(iconstr).eq.8) then
        rtot=0.d0
        do i=1,ndists(iconstr)
        at1=atmsconstr(iconstr,i)
        at2=atmsconstr(iconstr,i+1)

        rtot=dist(rclas(1,at1),rclas(2,at1),rclas(3,at1),
     .  rclas(1,at2),rclas(2,at2),rclas(3,at2))

        rtot=rtot+coef(iconstr,i)*rtot
c ndists
        enddo

        rt(iconstr)=rtot
        req=ro(iconstr)
        kf=kforce(iconstr)

        do i=1,ndists(iconstr)
        at1=atmsconstr(iconstr,i)
        at2=atmsconstr(iconstr,i+1)

c atom1: dr12
        dx=0.d0
        dy=0.d0
        dz=0.d0
          dx=(1.d0/rtot)*(rclas(1,at1)-rclas(1,at2))
          dx=2.d0*kf*(rtot-req)*dx
          dy=(1.d0/rtot)*(rclas(2,at1)-rclas(2,at2))
          dy=2.d0*kf*(rtot-req)*dy
          dz=(1.d0/rtot)*(rclas(3,at1)-rclas(3,at2))
          dz=2.d0*kf*(rtot-req)*dz
        fnew(1,1)=-dx*coef(iconstr,i)
        fnew(2,1)=-dy*coef(iconstr,i)
        fnew(3,1)=-dz*coef(iconstr,i)

c atom2: -fce over atom1
        fnew(1,2)=-fnew(1,1)
        fnew(2,2)=-fnew(2,1)
        fnew(3,2)=-fnew(3,1)

c adding fnew to fdummy
        fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)
        fdummy(1:3,at2)= fdummy(1:3,at2)+fnew(1:3,2)

c ndists
        enddo
c constrype
      endif

c writes variables for first step only
        if(iconstr.eq.1) then
        if(istep.eq.1) then
        write(*,'(/,A)')     'constr opt: variable constraint'
        write(*,'(A,i4)')    'icosntr  :', iconstr
        if(istepconstr.eq.1) then
        write(*,'(A,i4)')    'type     :', typeconstr(iconstr)
        if     (typeconstr(iconstr).eq.1) then
        write(*,'(A,4i4)')   'atoms    :', (atmsconstr(iconstr,i),i=1,4)
        elseif (typeconstr(iconstr).eq.2) then
        write(*,'(A,2i4)')   'atoms    :', (atmsconstr(iconstr,i),i=1,2)
        elseif (typeconstr(iconstr).eq.3) then
        write(*,'(A,3i4)')   'atoms    :', (atmsconstr(iconstr,i),i=1,3)
        elseif (typeconstr(iconstr).eq.4) then
        write(*,'(A,4i4)')   'atoms    :', (atmsconstr(iconstr,i),i=1,4)
        elseif (typeconstr(iconstr).eq.5) then
        write(*,'(A,4i4)')   'atoms    :', (atmsconstr(iconstr,i),i=1,4)
        elseif (typeconstr(iconstr).eq.6) then
        write(*,'(A,8i4)')   'atoms    :', (atmsconstr(iconstr,i),i=1,8)
        elseif (typeconstr(iconstr).eq.7) then
        write(*,'(A,5i4)')   'atoms    :', (atmsconstr(iconstr,i),i=1,5)
        elseif (typeconstr(iconstr).eq.8) then
        write(*,'(A,i4)')    'ndists   :', ndists(iconstr)
        do i=1,ndists(iconstr)
        write(*,'(A,F5.2)')  'coefs    :', coef(iconstr,i)
        enddo
        do i=1,ndists(iconstr)*2
        write(*,'(A,i4)')    'atoms    :', atmsconstr(iconstr,i)
        enddo
        endif
        write(*,'(A,F8.3)')  'constant :', kforce(iconstr)
        write(*,'(A,i4)')    'steps    :', nstepconstr
        write(*,'(A,F8.3)')  'initial  :', rini
        write(*,'(A,F8.3)')  'final    :', rfin
        endif
        write(*,'(A,F8.3)')  'targeted :', ro(iconstr)
        write(*,'(A,F8.3)')  'actual   :', rt(iconstr)
        endif
        else
        if(istep.eq.1) then
        write(*,'(/,A)')     'constr opt: fixed constraints'  
        write(*,'(A,i4)')    'icosntr  :', iconstr
        if(istepconstr.eq.1) then
        write(*,'(A,i4)')    'type     :', typeconstr(iconstr)
        if     (typeconstr(iconstr).eq.1) then
        write(*,'(A,4i4)')   'atoms    :', (atmsconstr(iconstr,i),i=1,4)
        elseif (typeconstr(iconstr).eq.2) then
        write(*,'(A,2i4)')   'atoms    :', (atmsconstr(iconstr,i),i=1,2)
        elseif (typeconstr(iconstr).eq.3) then
        write(*,'(A,3i4)')   'atoms    :', (atmsconstr(iconstr,i),i=1,3)
        elseif (typeconstr(iconstr).eq.4) then
        write(*,'(A,4i4)')   'atoms    :', (atmsconstr(iconstr,i),i=1,4)
        elseif (typeconstr(iconstr).eq.5) then
        write(*,'(A,4i4)')   'atoms    :', (atmsconstr(iconstr,i),i=1,4)
        elseif (typeconstr(iconstr).eq.6) then
        write(*,'(A,8i4)')   'atoms    :', (atmsconstr(iconstr,i),i=1,8)
        elseif (typeconstr(iconstr).eq.7) then
        write(*,'(A,5i4)')   'atoms    :', (atmsconstr(iconstr,i),i=1,5)
        elseif (typeconstr(iconstr).eq.8) then
        write(*,'(A,i4)')    'ndists   :', ndists(iconstr)
        do i=1,ndists(iconstr)
        write(*,'(A,F5.2)')  'coefs    :', coef(iconstr,i)
        enddo
        do i=1,ndists(iconstr)*2
        write(*,'(A,i4)')    'atoms    :', atmsconstr(iconstr,i)
        enddo
        endif
        write(*,'(A,F8.3)')  'constant :', kforce(iconstr)
        endif
        write(*,'(A,F8.3)')  'targeted :', ro(iconstr)
        write(*,'(A,F8.3)')  'actual   :', rt(iconstr)
        endif
        endif

c nconstr
        enddo 

c change units 
        rclas(1:3,1:natot)=rclas(1:3,1:natot)/0.529177d0
        end

c**************************************************************
      subroutine subconstr3(ro,rt,dr,E)

        implicit none
        double precision ro,rt,dr,E,eV,Ang
c        data eV /13.60580d0/
        data eV / 27.211396132d0/
        data Ang /0.529177d0/
        save eV,Ang 

	ro=ro+dr
	write(6,*)
	write(6,'(a)') 'constr opt: Optimized Reaction Coordinate'
	write(6,'(a,F10.5)') 'Reaction Coordinate (Ang) : ',rt
	write(6,'(a,F12.5)') 'Total Energy (eV) : ',E*eV
        end
c**************************************************************

