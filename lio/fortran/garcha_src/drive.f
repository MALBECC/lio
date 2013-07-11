c----------------------------------------------------------------
c READS AND DRIVES SUBROUTINE
c reads and prepare basis sets and geometry
c to perform electronic structure calculations
c using density functional theory and gaussian basis sets for
c the expansion of MO, charge density and exchange correlation
c potential
c Started 13 March 1992
c Up to now, it reads everything and builds up the basis
c repeats automatically 3p, 6d, 10d etc
c so is not necessary to do anything with that and
c input file is simple
c----------------------------------------------------------------
c
      subroutine drive(ng2,ngDyn,ngdDyn,P,scratch,vec)
       use latom 
c
      implicit real*8 (a-h,o-z)
      INCLUDE 'param'
      integer map(ntq)
      parameter (nng=500)
c
      character*65 title
      character*32 name1,whatis,stdbas
      character*4 date
      character*15 name2,name3,name4,name5,name6,solv,solv2
      character*4 ctype
      logical*4 exists,MEMO
      logical done(ntq),used,NORM,OPEN,ATRHO,DIRECT,VCINP,SHFT
      logical TMP,TMP2,dens,EXTR,write,SVD,ANG,field,field1
      logical Ex,Coul,Scf1,Prop,GRAD,BSSE,integ,SVD1,sol,tipe
      logical exter,exter1,resp1,popf,parsearch
c
      logical OPEN1
      logical dens1,integ1,sol1,free,free1
c
      dimension OCC(40),oc2(400),ATCOEF(100*ng0),ighost(ntq),
     > ighost1(ntq)
      dimension ncf(nng),lt(nng)
      namelist /int0/ title,natom,nopt,NORM,ANG,Ex,Coul,Scf1,Prop,
     >                field,sol,resp1
      namelist /intMD/ h,nsteps,istart
      namelist /scfinp/ OPEN,NMAX,NCO,Nunp,ATRHO,VCINP,DIRECT,
     >     EXTR,SHFT,SHI,IDAMP,GOLD,told,write,MEMO,rmax,rmaxs,parsearch
      namelist /RHOINP/ NAO,NGF,OCC,ATCOEF
      namelist /EXCH/ Iexch,integ,dens,igrid,igrid2
      namelist /COULx/ SVD,iconst
      namelist /GEO/ FTOL,imin,ibrent,delta,inmod,thermo,TEMP,sigma
      namelist /Propx/ idip,ipop,ispin,icharge,map,popf
      namelist /BSSE1/ nco1,nunp1,open1,ighost1
      namelist /fieldx/ a0,epsilon,exter,Fx,Fy,Fz
      namelist /solx/ Nsol,natsol,solv,free
      namelist /rfield/ F
      common /intg1/ e_(50,3),wang(50),Nr(0:54)
      common /intg2/ e_2(116,3),wang2(116),Nr2(0:54),e3(194,3),
     > wang3(194)
c
      dimension x(nt),y(nt),z(nt),Pm(nt),vx(nt),vy(nt),vz(nt)
c
      dimension Iz(nt),r(nt,3),nnat(nt),v(nt,3),isotop(54)
      equivalence (r(1,1),x(1)),(r(1,2),y(1)),(r(1,3),z(1))
      equivalence (v(1,1),vx(1)),(v(1,2),vy(1)),(v(1,3),vz(1))
c
      dimension Em(ntq+nss),Rm(ntq+nss),pc(nss)
      dimension alpha(nss)
c
      COMMON /TABLE/ STR(880,0:21)
      common /Sys/ SVD1,iconst1
      common /dyn/ h1,nsteps1,Pm,Nscale
      common/Ngeom/ ngeo
      common /ENum/ GRAD
      common /bsse/ bsse,ighost
      common /propt/ idip1,ipop1,ispin1,icharge1,map1(ntq)
      common /masses/ xmass(216)
c
      common /sol1/ Nsol1,natsol1,alpha,Em,Rm,pc,sol1,free1
      common /radii/ Rm2(0:54)
c
c Everything is dimensioned for 2 basis, normal and density
c ncf, lt,at,ct parameters for atomic basis sets
      dimension at(nng),ct(nng),nshell(0:4)
      dimension Num(0:3),nlb(ng),nld(ngd),nshelld(0:4)
 
c
c auxiliar , and debuggings
c ------------------
c parameters for 2 basis sets, normal, density 
c ng maximum # of contracted functions , nl maximum # of primitives in
c a contraction
c c, cd ce , linear coefficients
c a , ad, ae exponents
c Nuc , indicates center , ncont # of primitives in the correspondent
c contraction, and nx ny nz exponents of x , y and z in the cartesian
c gaussian
c ------------------
      dimension c(ng,nl),a(ng,nl),Nuc(ng),ncont(ng)
      dimension cx(ngd,nl),ax(ngd,nl),Nucx(ngd),ncontx(ngd)
      dimension cd(ngd,nl),ad(ngd,nl),Nucd(ngd),ncontd(ngd)

      integer ii(ng),iid(ngd)
c
c Dimensions for Overlap and Fock matrices
c
      dimension P(ng2),scratch(ngDyn,3*ngDyn),vec(ngdDyn,ngdDyn)
c
      parameter (pi=3.14159265358979312D0)
c Angular momenta : up to f functions ( could be easily extended if
c necessary)
c
c     common /fit/ Nang,dens,integ,Iexch,igrid,igrid2
      common /fit/ Nang,dens1,integ1,Iexch1,igridc,igrid2c
      common /cav/ a01,epsilon1,field1,exter1,Fx1,Fy1,Fz1
c
      common /index/ ii,iid
c
      Data Num /1,3,6,10/
c
c --- Defaults ---------------------------------------------------------
c
c NORM true , expansion in normalized gaussians, so normalization factor
c included in coefficients
c default most stable isotopes, for mass asignation
      do i = 1,54
        isotop(i) = 1
      enddo
c
      Nunp=0
      iconst=0
      idip=1
      ispin=1
      ipop=1
C
	icharge=0
C
      NMAX=50
      write=.false.
      NORM=.true.
      SVD=.false.
      title='               '
      natom=1
      nopt=0
      Nang=50
      Iexch=2
      dens=.false.
      OPEN=.false.
      DIRECT=.true.
      ATRHO=.true.
      GOLD=1.D0
      SHFT=.false.
      ANG=.false.
      SHI=1.D0
      IDAMP=0
      told=1.0D-4
      EXTR=.false.
      istart=0
      SVD=.false.
      Ex=.true.
      Coul=.false.
      Scf1=.true.
      Prop=.false.
      igrid=1
      igrid2=0
      iforce=0
      ngeo=0
      FTOL=1.0D-04
      delta = 0.01D0
      inmod = 0
      imin=0
      ibrent=0
      GRAD=.true.
      BSSE=.false.
      integ=.false.
      scale=1.0D0
      scale2=1.0D0
      thermo=0
      TEMP=298.0D0
      sigma=1
      a0=1000.0
      epsilon=1.0D0
      field=.false.
      popf=.false.
c
      h=0.15D0
      nsteps=10
      Nscale=20000
c
      sol=.false.
      free=.false.
      rmax=16.0D0
      rmaxs=3.0D0

c
      parsearch=.false.

      do n=1,ntq
       map(n)=n
      enddo
c calls generator of table for incomplete gamma functions
c
       call GENERF
       call GENERFS
       call GRID
c-----------------------------------------------------------------------
c reads input file
      write(*,*) ' Name of the input file ?'
      read(*,100) name1
c
c
      inquire(file=name1,exist=exists)
      if (.not.exists) then
      write(*,*) 'ERROR CANNOT FIND INPUT FILE ON UNIT 1'
      pause
      else
c  name of output file
      ikk=1
      do 1 while (name1(ikk:ikk).ne.' ')
       ikk=ikk+1
 1    continue
c
      ikk=ikk-1
      name2=name1(1:ikk)//'.out'
c
      open(unit=1,file=name1,iostat=ios)
      open(unit=2,file=name2)
      endif
c-------------------------------------------------------
      date='date'
      write(*,*) 'JOB STARTED '
      call system(date)
c
c --- Reads input file
c
c nopt 0 for SCF , 1 for MD-SA , 2 for Standard Geometry Optimization
c nopt 3 for SCF-nuclear MD
c natom , all in namelist int0
c
      read(1,nml=int0)
      write(2,nml=int0)
c
      nopt1=nopt
c
      if (nopt.eq.0) then
      name3=name1(1:ikk)//'.orb'
      open (unit=29,file=name3)
      endif
c
      sol1=sol
      if (sol) then
       read(1,nml=solx)
       write(2,nml=solx)
       free1=free
      endif
c
c-------------------------------------------------------
c opens file with solvent parameters and reads them
      if (sol) then
       ik2=1
      do  while (solv(ik2:ik2).ne.' ')
       ik2=ik2+1
      enddo
c
      ik2=ik2-1
      solv2=solv(1:ik2)//'.param'
c
      open(unit=4,file=solv2)
c
      do n=1,natsol
       read(4,*) Em(natom+n),Rm(natom+n),pc(n),alpha(n)
      enddo
c ----
c now it should read also solute Lennard Jones parameters
c it could be done automatically in the future
c
      name6=name1(1:ikk)//'.param'
      open(unit=7,file=name6)
c
      do n=1,natom
       read(7,*) Em(n),Rm(n)
      enddo
c
      endif
c-------------------------------------------------------
      Nsol1=Nsol
      natsol1=natsol
c
      if (Prop) then
      read(1,nml=propx)
      write(2,nml=propx)
      endif
c
      if (nopt.eq.4) then
       read(1,nml=BSSE1)
       do n=1,natom
       ighost(n)=ighost1(n)
       enddo
      endif
c
c geometry, in the MD case also initial velocities
c  and nuclear masses needed
c
c initial coordinates
      do 10 i=1,natom
       read(1,*) Iz(i),x(i),y(i),z(i)
c
       indmass = (Iz(i)-1)*4 + isotop(Iz(i))
       Pm(i) = xmass(indmass)
c -----------------------------------------------------
       done(i)=.false.
       write(2,500) Iz(i),x(i),y(i),z(i)
  10  continue
c
c only in case of mixed calculations
c
      if (sol) then
       ntatom=natom+nsol*natsol
c
       do 19 i=natom+1,ntatom
       read(1,*) Iz(i),x(i),y(i),z(i)
c
       indmass = (Iz(i)-1)*4 + isotop(Iz(i))
       Pm(i) = xmass(indmass)
       write(2,501) Iz(i),x(i),y(i),z(i)
 19    continue
c
      else
       ntatom=natom
      endif
c
c
       if (ANG) then
         do 17 i=1,ntatom
          x(i)=x(i)/0.529177D0
          y(i)=y(i)/0.529177D0
          z(i)=z(i)/0.529177D0
  17     continue
       endif
c----------------
c REACTION FIELD CASE
        if (field) then
        read(1,nml=fieldx)
        write(2,nml=fieldx)
        endif
c ELECTRICAL RESPONSE CASE
        if (resp1) then
         read(1,nml=rfield)
         write(2,nml=rfield)
        endif
c----------
c

       if (nopt.eq.1.or.nopt.eq.3) then
c time step , number of steps and masses
       read(1,nml=intMD)
       if (nopt.eq.3) then
        istart0=istart
        istart=1
       endif
       write(2,nml=intMD)
       endif
c
       h1=h
       nsteps1=nsteps
c
       if (nopt.eq.2) then
        read(1,nml=GEO)
        write(2,nml=GEO)
        if(inmod.ne.0) then
          name3=name1(1:ikk)//'.nmod'
          name4=name1(1:ikk)//'.freq'
          name5=name1(1:ikk)//'.vib'
          open(unit=45,file=name3)
          open(unit=46,file=name4)
          open(unit=44,file=name5)
        endif
       endif
c
c -------------------------------------------------------------
c for each kind of atom, given by Z
c
c Basis set for MO expansion
c
c reads whatis: Gaussian for the time being, 
c      later stdbas will be included
c reads iatom, atomic number , nraw number of primitive gaussians
c , counting 3 p 6 d .., and ncon number of contractions counted in the
c same way
c if coefficients correspond to normalized gaussians
c normalization factor will  multiply coefficients, so in  the program
c expressions for non-normalized gaussians are used.
c an option is included for using non-normalized gaussians
c !The contractions don't need to be normalized, it will be done
c automatically
c
c
      read(1,100) whatis
      write(2,100) whatis
c
c
      No=0
      Nd=0
      Ne=0
      NBAS=0
    	M=0
	    Mdd=0

        allocate(natomc(natom),nnps(natom),nnpp(natom),nnp(natom))
        allocate(nnpd(natom),nns(natom),nnd(natom),atmin(natom))
        allocate(jatc(natom,natom))
       
         do i=1,natom
         natomc(i)=0 
          enddo
c
c  BASIS SETS -------------------------------------------------
      do 25 while (whatis.ne.'endbasis')
c
      NBAS=NBAS+1
c signals if a basis set was not used
      used=.false.
      atmint=100000.
      read(1,*) iatom,nraw,ncon
      write(2,600) iatom,nraw,ncon
c
c reads contraction scheme. The value for p,d ,f should not be repeated
c 3 ,6 , 10 .....   times
c reads also angular momentum for each of the contractions
c 0 for s 1 for p, etc
c
       
      read(1,*) (ncf(i),i=1,ncon)
      write(2,*) (ncf(i),i=1,ncon)
      read(1,*) (lt(i),i=1,ncon)
      write(2,*) (lt(i),i=1,ncon)
c
c loop over all primitives, no repeating p, d
      do 30 i=1,nraw
       read(1,*) at(i),ct(i)
       write(2,700) at(i),ct(i)
c       write(*,*) atmint,at(i)
        if(at(i).lt.atmint) atmint=at(i)

 30   continue



c
      do 35 j=1,natom
       if ((Iz(j).eq.iatom).and.(.not.done(j))) then
       nnat(NBAS)=nnat(NBAS)+1
       done(j)=.true.
       used=.true.
       atmin(j)=atmint

c  cosas que puso nano para "atomos cerca"

        nns(j)=0
        nnp(j)=0
        nnd(j)=0        

       do kkk=1,ncon
        if (lt(kkk).eq.0) nns(j)=nns(j)+Num(lt(kkk))
        if (lt(kkk).eq.1) nnp(j)=nnp(j)+Num(lt(kkk))
        if (lt(kkk).eq.2) nnd(j)=nnd(j)+Num(lt(kkk))
       enddo

c      write(*,*) 'nns y etc',nns(j),nnp(j),nnd(j),nnps(j)
c     > ,nnpp(j),nnpd(j)

c  M stores # of contractions
       index=0
       do 36 k=1,ncon
c
        M=M+Num(lt(k))
c nshell gives the # of functions s, p, d  etc
        nshell(lt(k))=nshell(lt(k))+Num(lt(k))
c
       do 36 l2=1,Num(lt(k))
c
       No=No+1
c
c normalization
c
        if (NORM) then
        do 37 l=1,ncf(k)
c 
        index=index+1
        goto (70,80,90) lt(k)+1
c
 70     xnorm=sqrt((2.D0*at(index)/pi)**3)
        xnorm=sqrt(xnorm)   
        c(No,l)=ct(index)*xnorm
        a(No,l)=at(index)
      goto 37
c
 80     xnorm=sqrt((2.D0*at(index)/pi)**3)*4.D0*at(index)
        xnorm=sqrt(xnorm)
        c(No,l)=ct(index)*xnorm
        a(No,l)=at(index)
      goto 37
c
 90     xnorm=sqrt((2.D0*at(index)/pi)**3)*(4.D0*at(index))**2
        xnorm=sqrt(xnorm)
        c(No,l)=ct(index)*xnorm
        a(No,l)=at(index)
      goto 37
c
 37   continue
c
      else
c no normalization case
c
       do 38 l=1,ncf(k)
        index=index+1
c
        c(No,l)=ct(index)
        a(No,l)=at(index)
 38   continue
c
      endif
c repeat the index only for p,d and f
       if (l2.ne.Num(lt(k))) then
        index=index-ncf(k)
       endif
c
      Nuc(No)=j
      ncont(No)=ncf(k)
      nlb(No)=lt(k)
 36   continue
      endif
 35   continue
c     
c
      if (.not.used) then
       write(*,200) iatom
      endif
c
c Exactly the same should be repeated for charge density
c and exchange correlation potential basis sets
c
c CHARGE DENSITY --------------------------------------------------
c
      read(1,*) iatom,nraw,ncon
      write(2,*) iatom,nraw,ncon
c
c reads contraction scheme. The value for p,d ,f should not be repeated
c 3 ,6 , 10 .....   times. Reads also angular type , 
c 0 for s , 1 for p etc
      read(1,*) (ncf(i),i=1,ncon)
      write(2,*) (ncf(i),i=1,ncon)
      read(1,*) (lt(i),i=1,ncon)
      write(2,*) (lt(i),i=1,ncon)
c     
c
c loop over all primitives, repeating p, d
      do 40 i=1,nraw
       read(1,*) at(i),ct(i)
       write(2,700) at(i),ct(i)
 40   continue
c
      do 45 j=1,natom
       if (Iz(j).eq.iatom) then
c
c Mdd stores # of contractions in final basis, counting all possibilities
c for p , d etc
c
       index=0
       do 46 k=1,ncon
c
       Mdd=Mdd+Num(lt(k))
 
       nshelld(lt(k))=nshelld(lt(k))+Num(lt(k))
c
       do 46 l2=1,Num(lt(k))
c
       Nd=Nd+1
c
        if (NORM) then
       do 47 l=1,ncf(k)
        index=index+1
c
        goto (71,81,91) lt(k)+1
c
 71     xnorm=sqrt((2.D0*at(index)/pi)**3)
        xnorm=sqrt(xnorm)
        cd(Nd,l)=ct(index)*xnorm
        ad(Nd,l)=at(index)
c       ad(Nd,l)=2.D0*at(index)
      goto 47
c
 81     xnorm=sqrt((2.D0*at(index)/pi)**3)*4.D0*at(index)
        xnorm=sqrt(xnorm)
        cd(Nd,l)=ct(index)*xnorm
c       ad(Nd,l)=2.D0*at(index)
        ad(Nd,l)=at(index)
      goto 47
c
 91     xnorm=sqrt((2.D0*at(index)/pi)**3)*(4.D0*at(index))**2
        xnorm=sqrt(xnorm)
        cd(Nd,l)=ct(index)*xnorm
c       ad(Nd,l)=2.D0*at(index)
        ad(Nd,l)=at(index)
      goto 47
c
 47   continue
c
      else
c
c no normalization case
c
       do 48 l=1,ncf(k)
        index=index+1
c
        cd(Nd,l)=ct(index)
c       ad(Nd,l)=2.D0*at(index)
        ad(Nd,l)=at(index)
 48   continue
c
      endif
c
c repeat the index only for p,d and f, criterium l2<max(l2)
       if (l2.ne.Num(lt(k))) then
        index=index-ncf(k)
       endif
c
      Nucd(Nd)=j
      ncontd(Nd)=ncf(k)
      nld(Nd)=lt(k)
c
 46   continue
      endif
 45   continue
c
      read(1,100) whatis
      write(2,100) whatis
 25   end do
c----- DIMENSION CONTROLS ------------------------------------
c
c       do i=1,natom
c      write(*,*) 'atmin',i,atmin(i)
c      enddo


      iprob=0
      if (M.gt.ngDyn.or.M.gt.ng) then
       write(*,*) 'DIMENSION PROBLEMS WITH BASIS SET SIZE PARAMETER NG'
       write(*,*) 'NUMBER BASIS FUNCTIONS =',M,'<',ngDyn
       iprob=1
      endif
c
      if (Mdd.gt.ngdDyn.or.Mdd.gt.ngd) then
       write(*,*) 'DIMENSION PROBLEMS WITH AUXILIAR BASIS PARAMETER NGD'
       write(*,*) 'NUMBER AUXILIAR BASIS FUNCTIONS =',Mdd,'<',ngdDyn
       iprob=1
      endif
c
      if (natom.gt.nt) then
       write(*,*) 'DIMENSION PROBLEMS WITH NUMBER OF ATOMS PARAMETER NT'
       write(*,*) 'NUMBER OF ATOMS =',natom,'<',nt
       iprob=1
      endif
c
c -------------------------------------------------------------
c -- Initial guess -------------------------------------------
c
c
c--- changes basis order to one in which the order is given by
c the angular type :all s first, then all p, then all d ,......
c also converts into normalized gaussians, including normalization
c factor into linear coefficients
c for d shell, in order to keep the same contraction for all the shell
c x^2, y^2 and z^2 are normalized to 3, later in the integral part,
c this should be considered
c
      is=1
      ip=nshell(0)+1
      id=nshell(0)+nshell(1)+1
c
c standard basis set  ---------------------------------------------
c
c loop over the total # of basis
      do 210 i=1,M
c
      l=nlb(i)+1
      goto (11,22,33) l
c
c s case
  11  continue
c   
      Nucx(is)=Nuc(i)
      ii(is)=i
      ncontx(is)=ncont(i)
c  
      do 12 j=1,ncontx(is)
c
       cx(is,j)=c(i,j)
  12   ax(is,j)=a(i,j)
c
      is=is+1
      goto 44   
c
c p case
  22  continue
c   
      Nucx(ip)=Nuc(i)
      ii(ip)=i
      ncontx(ip)=ncont(i)
c  
      do 13 j=1,ncontx(ip)
c
       cx(ip,j)=c(i,j)
  13   ax(ip,j)=a(i,j)
c
      ip=ip+1
      goto 44
c
c d case
  33  continue
c   
      Nucx(id)=Nuc(i)
      ii(id)=i


      ncontx(id)=ncont(i)
c  
      do 14 j=1,ncontx(id)
c
       cx(id,j)=c(i,j)
  14   ax(id,j)=a(i,j)
c
      id=id+1
      goto 44
c
  44  continue
 210  continue
c
c final normalization for d
c 3 factor for x^2, y^2 and z^2
c is in the integral part, so shell structure could be used
c
c
c Now, it has to put temporary things into the real one
c
      do 300 i=1,M
       Nuc(i)=Nucx(i)
       ncont(i)=ncontx(i)
c
       do 300 j=1,ncont(i)
        c(i,j)=cx(i,j)
        a(i,j)=ax(i,j)
 300   continue
c
        
c same process, but for electronic density basis set ------------
c
      is=1
      ip=nshelld(0)+1
      id=nshelld(0)+nshelld(1)+1
c
c loop over the total # of basis
      do 310 i=1,Mdd
c
      l=nld(i)+1
      goto (111,222,333) l
c
c s case
 111  continue
c   
      Nucx(is)=Nucd(i)
      iid(is)=i
      ncontx(is)=ncontd(i)
c  
      do 128 j=1,ncontx(is)
c
c
       cx(is,j)=cd(i,j)
  128  ax(is,j)=ad(i,j)
c
      is=is+1
      goto 444   
c
c p case
 222  continue
c   
      Nucx(ip)=Nucd(i)
      iid(ip)=i
      ncontx(ip)=ncontd(i)
c  
c
      do 213 j=1,ncontx(ip)
c
       cx(ip,j)=cd(i,j)
 213   ax(ip,j)=ad(i,j)
c
      ip=ip+1
      goto 444
c
c d case
 333  continue
c   
      Nucx(id)=Nucd(i)
      iid(id)=i
      ncontx(id)=ncontd(i)
c  
      do 145 j=1,ncontx(id)
c
c
       cx(id,j)=cd(i,j)
  145  ax(id,j)=ad(i,j)
c
      id=id+1
      goto 444
c
 444  continue
 310  continue
c
c final normalization for d
c 3 factor for x^2, y^2 and z^2
c is in the integral part, so shell structure could be used
c
c
c Now, it has to put temporary things into the real one
c
      do 301 i=1,Mdd
       Nucd(i)=Nucx(i)
       ncontd(i)=ncontx(i)
c
       do 301 j=1,ncontd(i)
        cd(i,j)=cx(i,j)
        ad(i,j)=ax(i,j)
 301   continue
c
c------------------------------------------------------------------
c POINTERS --------------------------------------------
c
      MM=M*(M+1)/2
      MMd=Mdd*(Mdd+1)/2
c first P
      M1=1
c now F alpha
      M3=M1+MM
c now S, F beta also uses the same position after S was used
      M5=M3+MM
c now G
      M7=M5+MM
c now Gm
      M9=M7+MMd
c now H
      M11=M9+MMd
c W ( eigenvalues ), also this space is used in least squares
      M13=M11+MM
c aux ( vector for ESSl)
      M15=M13+M
c Least squares
      M17=M15+MM
c vectors of MO alpha
      M18=M17+MMd
c vectors of MO beta
c------------------------------------------------------------
c
c
c Density matrix  construction - For closed shell only
c
c
      if (Scf1) then
      read(1,nml=scfinp)
      endif

c
c variables defined in namelist cannot be in common ?
      OPEN1=OPEN
      NCO2=NCO
      Nunp2=Nunp
c
       idip1=idip
       ipop1=ipop
       icharge1=icharge
       ispin1=ispin
       do n=1,ntq
        map1(n)=map(n)
       enddo
c
       if (ANG) then
       a01=a0/.529177D0
       else
       a01=a0
       endif
       
       
       epsilon1=epsilon
       field1=field
       exter1=exter
       Fx1=Fx
       Fy1=Fy
       Fz1=Fz
c
      if (Ex) then
      read(1,nml=EXCH)
      endif
c
      dens1=dens
      integ1=integ
      Iexch1=Iexch
      igridc=igrid
      igrid2c=igrid2
c
      if ((Iexch.ge.4).and.(.not.(integ)).and.(.not.(dens))) then
       write(*,*) 'OPTION SELECTED NOT AVAILABLE'
c      pause
      endif
c
      if ((Iexch.eq.2).and.(OPEN)) then
       write(*,*) 'OPTION SELECTED NOT AVAILABLE YET'
c      pause
      endif
c
      if (Coul) then
      read(1,nml=COULx)
      endif
c
      SVD1=SVD
      iconst1=iconst
c DIMENSION TESTS -----------------------------------------------
c
      Ndim=5*M*(M+1)/2+3*Mdd*(Mdd+1)/2+M+M*NCO+M*Ngrid
      if (Ndim.gt.ng2) then
       write(*,*) 'DIMENSION PROBLEMS WITH DYNAMICAL VECTOR NG2'
       iprob=1
      endif
c
      if (iprob.eq.1) then
       pause
      endif
c
c -------------------------------------------------------------
c case for initial guess given in input -------------------------
      if (VCINP) then
c closed shell
c     
      if (.not.OPEN) then
c reads vectors of MO coefficients, basis is in the same order as 
c the given in input
      do 120 l=1,M
 120   read(1,*) (vec(l,n),n=1,NCO)
c
c puts vectors in dynamical allocation (to be used later)
c
       kk=0
      do 121 k=1,NCO
      do 121 i=1,M
       kk=kk+1
 121   P(M18+kk-1)=vec(ii(i),k)

      do 130 i=1,M
      do 130 j=1,M
c
      do 130 k=1,NCO
       scratch(i,j)=scratch(i,j)+2.0D0*vec(i,k)*vec(j,k)
 130  continue
c
c open shell case
      else
c
      NCOa=NCO
      NCOb=NCO+Nunp
      M18b=M18+M*NCOa
c alpha
      do 220 l=1,M
 220    read(1,*) (vec(l,n),n=1,NCOa)
c

       kk=M18-1
      do 221 k=1,NCOa
      do 221 i=1,M
       kk=kk+1
 221   P(kk)=vec(ii(i),k)
c
c Density Matrix
c
      do 331 i=1,M
      do 331 j=1,M
      scratch(i,j)=0.0D0
c
      do 139 k=1,NCOa
       scratch(i,j)=scratch(i,j)+vec(i,k)*vec(j,k)
 139  continue
 331  continue
c
c beta
      do 229 l=1,M
 229    read(1,*) (vec(l,n),n=1,NCOb)
c
       kk=M18b-1
      do 223 k=1,NCOb
      do 223 i=1,M
       kk=kk+1
 223   P(kk)=vec(ii(i),k)
c 
c Density Matrix
c
      do 335 i=1,M
      do 335 j=1,M
c
      do 239 k=1,NCOb
       scratch(i,j)=scratch(i,j)+vec(i,k)*vec(j,k)
 239  continue
 335  continue
c
      endif
      endif
c
c density matrix kept temporarily in S
c end of case in which density matrix is explicitly given in
c input
c-----------------------------------------------------------------
c
c case for initial guess constructed from atomic densities -------
      if (ATRHO) then
c
      k1=0
      l1=0
      NN=0
      do 152 i=1,NBAS
c
        read(1,nml=RHOINP)
c
      do 152 j=1,nnat(i)
c
c
      do 150 k=1,NAO
       k1=k1+1
c
      do 151 l=1,NGF
       l1=NN+l
c
       kl=(k-1)*NGF+l
c
       vec(l1,k1)=atcoef(kl)
       oc2(k1)=OCC(k)
 151  continue
c
c
 150  continue
c
      NN=NN+NGF

 152  continue
c
c 
c------------------------------------------------------------
c
c
c S used as scratch array here
      do 118 i=1,M
      do 118 j=1,M
c
      scratch(i,j)=0.
      do 109 l=1,k1
  109 scratch(i,j)=scratch(i,j)+oc2(l)*vec(i,l)*vec(j,l)
  118 continue
c
c approximate vectors construction
c
      if (.not.OPEN) then
       kk=M18-1
      do 430 k=1,NCO
      do 430 i=1,M
       kk=kk+1
 430  P(kk)=vec(i,k) 
c
      else
c
      NCOa=NCO
      NCOb=NCO+Nunp
      M18b=M18+M*NCOa
c
      kk=M18-1
      do 431 k=1,NCOa
      do 431 i=1,M
       kk=kk+1
 431  P(kk)=vec(ii(i),k)
c
      kk=M18b-1
      do 432 k=1,NCOb
      do 432 i=1,M
       kk=kk+1
 432  P(kk)=vec(ii(i),k)
      endif
c

      endif
c density matrix stored temporarily in S, then it should be changed
c according to the shell ordering of the basis set and kept in P
c------ end of option atomic densities -----------------
c
c
c changes to the shell order ( s , p, d....)
      k=0
      do 119 j=1,M
      do 119 i=j,M
       k=k+1
       P(k)=Scratch(ii(i),ii(j))
 119  continue
c
      k=0
      do 127 j=1,M
      do 127 i=j,M
       k=k+1
       if (i.ne.j) then
       P(k)=P(k)*2.D0
       endif
 127   continue
c
c---- reads exchange fit data -------------
c
      TMP=ATRHO
      TMP2=VCINP
      ATRHO=.FALSE.
      VCINP=.TRUE.
c this only for NH4Cl--------
c----------------------------
      if (Scf1) then
      write(2,nml=scfinp)
      endif
c
      if (ex) then
      write(2,nml=EXCH)
      endif
c
      if (coul) then
      write(2,nml=COULx)
      endif
c
      ATRHO=TMP
      VCINP=TMP2
c--------------------------
c-------------------------------------------
c reads geometry and velocities in case of MD
c
c in case of SCF calculation, geometry of restart file should be same
c this is only for MD
      if ((nopt.eq.3).and.(istart0.eq.1)) then
      do ki=1,ntatom
       read(1,*) Iz(ki),x(ki),y(ki),z(ki)
      enddo
c
      do  kk=1,ntatom
       read(1,*) vx(kk),vy(kk),vz(kk)
       if (ANG) then
          tmpr=1.0D0/0.529177D0
          x(kk)=x(kk)*tmpr
          y(kk)=y(kk)*tmpr
          z(kk)=z(kk)*tmpr
          vx(kk)=vx(kk)*tmpr
          vy(kk)=vy(kk)*tmpr
          vz(kk)=vz(kk)*tmpr
      endif
      enddo
c
      endif
c
      M2=2*M
c
c Reaction Field - Calculation of Molecular Volume and Size of cavity.
c algorithm given in : JACS 114 1645 (1992).
      if (field) then
c
       if (exter) then
       write (*,*) 'CALCULATION WITH AN EXTERNAL FIELD'
c
       else
       if (a0.ge.1000.0D0) then
       call vol(NORM,natom,r,Nuc,Iz,M,Mdd,ncont,nshell,
     >            c,a,M18,NCO,NCO+Nunp,OPEN1,P,a01)
      endif
      write (*,*) 'REACTION FIELD CALCULATION WITH'
      write (*,320) a01,epsilon1
      write(*,*)

      endif
c
      endif
c
       if (resp1) then
        nopt=8
        exter1=.true.
        field1=.true.
        write(*,650) F
        call resp(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,
     >  nshell,c,a,Nucd,Mdd,ncontd,nshelld,cd,ad,P,scratch,vec,F,
     >   nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
        return
       endif
c calculation of fragment orbital populations
       if (popf) then
       call popu(NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >         Mdd,P,scratch,vec,ngdDyn,Nucx,
     >   nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
       stop
       endif
c-----------------------------
#ifdef G2G
c------- G2G Initialization ---------------------
      call g2g_parameter_init(NORM, natom, ntq, ng, r, Rm2, Iz, Nr, Nr2,
     > Nuc, M, ncont, nshell, c, a, P, M18, M5, NCO, nopt, Iexch,
     > e_, e_2, e3, wang, wang2, wang3)
c       write(99,*) NORM, natom, ntq, ng, r, Rm2, Iz, Nr, Nr2,
c     > Nuc, M, ncont, nshell, c, a, P, M18, M5, NCO, nopt, Iexch,
c     > e_, e_2, e3, wang, wang2, wang3
c      write(98,*) Rm2
#endif
c nopt 0 static SCF calculation --------------------------------------
      if (nopt.eq.0) then

#ifdef G2G
			write(*,*) 'primera carga de posiciones'
			call g2g_reload_atom_positions(igrid2)
#endif

      if (parsearch) then
       call g2g_solve_groups(1, E, 0)
       goto 101
      endif
     



      if (OPEN) then
c Nunp : number of unpaired electrons
c
      call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >         Nucd,Mdd,ncontd,nshelld,cd,ad,P,scratch,vec,E,
     >   nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
c
      else
      call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >         Nucd,Mdd,ncontd,nshelld,cd,ad,P,scratch,vec,E,
     >   nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)


c
c
      endif
      endif
c
c -- BSSE option -------------------------------------------------
      if (nopt.eq.4) then
c 
      write(*,*)
      write(*,*) ' SCF ENERGY CALCULATION WITH GHOST ATOMS (BSSE)'
c
c it is convenient to use as a starting guess , the vectors from
c a full system calculation. The terms corresponding to the ghost
c atoms is set to 0
c
c
      k=0
      do j=1,M
      do i=j,M
       k=k+1
       if (ighost(Nuc(i))*ighost(Nuc(j)).eq.0) then
       P(k)=0.0D0
       endif
      enddo
      enddo
c
      BSSE=.true.
      NCO2=nco1
      Nunp2=Nunp1

c

      if (OPEN1) then
       OPEN=.true.
      call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >         Nucd,Mdd,ncontd,nshelld,cd,ad,P,scratch,vec,E,
     >   nopt,OPEN,NMAX,NCO2,ATRHO,VCINP,SHFT,Nunp2,GOLD,told,write)
c
      else
      call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >         Nucd,Mdd,ncontd,nshelld,cd,ad,P,scratch,vec,E,
     >   nopt,OPEN,NMAX,NCO2,ATRHO,VCINP,SHFT,Nunp2,GOLD,told,write)

      endif
      endif
c
c------------------------------------------------------
c Standard (SCF) molecular dynamics option
      if (nopt.eq.3) then
      if (istart.eq.0) then
       do kk=1,ntatom
        vx(kk)=0.D0
        vy(kk)=0.D0
        vz(kk)=0.D0
       enddo
      endif
c
      call MD2(MEMO,NORM,natom,Iz,r,v,Nuc,M,ncont,nshell,c,a,
     >     Nucd,Mdd,ncontd,nshelld,cd,ad,P,scratch,vec,E,name1,ikk,
     >   nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
c
      endif
c
c--- Geometry optimization option
c
      if (nopt.eq.2) then
       call geom(MEMO,FTOL,imin,ibrent,NORM,natom,Iz,r,Nuc,M,ncont,
     >   nshell,c,a,Nucd,Mdd,ncontd,nshelld,cd,ad,P,scratch,vec,E,Pm,
     >           delta,inmod,name1,ikk,thermo,TEMP,sigma,
     >   nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
c
       endif
c
c-----------------------------------------------------------------------
c ---- Output results --------
c vector of molecular orbitals is going to be written in SCF subroutine
c
      date='date'
      write(*,*) 'JOB FINISHED'
      call system(date)

101   DEALLOCATE (natomc,nnps,nnpp,nnpd,nns)
      deallocate (nnd,nnp,atmin,jatc)

c
c---------------------------------------------------
c---------------------------------------------------
 100  format (A8)
 200  format ('basis set corresponding to Z ',I3,' was not used')
 400  format ('not implemented for open shell yet')
 500  format (i3,3x,F11.6,2x,F11.6,2x,F11.6)
 501  format (i3,3x,F11.6,2x,F11.6,2x,F11.6, ' CLASSICAL')
 600  format (3(i2,2x))
 650  format ('Electric response calculation F =',F7.4)
 700  format (F15.7,3x,F9.6)
 320  format (' Cavity Size (a.u)',F9.3,'     Dielectric Constant',F7.2)
c
      return
      end

