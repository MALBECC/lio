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
      subroutine drive(ng2,ngDyn,ngdDyn)
       use garcha_mod

c
      implicit real*8 (a-h,o-z)
      logical Exx, parsearch
      character*255 int_basis_file, fit_basis_file
      character*255 liohome
      character*255 inp_line
      character inp_char
c      namelist /scfinp/ OPEN,NMAX,Nunp,ATRHO,VCINP,DIRECT,
c     >EXTR,SHFT,SHI,IDAMP,GOLD,told,write1,MEMO,rmax,rmaxs,predcoef,
c     >idip,writexyz,intsoldouble,watermod,DIIS,ndiis,dgtrig
c      namelist /EXCH/ Iexch,integ,dens,igrid,igrid2
c      namelist /COULx/ SVD,iconst
c      namelist /GEO/ FTOL,imin,ibrent,delta,inmod,thermo,TEMP,sigma
c      namelist /Propx/ ipop,ispin,icharge,popf,parsearch!,map
c      namelist /BSSE1/ nco1,nunp1,open1,ighost1
c      namelist /fieldx/ a0,epsilon,exter,Fx,Fy,Fz
c      namelist /solx/ Nsol,natsol,solv,free
c      namelist /rfield/ F
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
c
c Dimensions for Overlap and Fock matrices
c
c
c Angular momenta : up to f functions ( could be easily extended if
c necessary)
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
      ngeo=0
      MEMO=.true.
      iconst=0
      idip=1
      ispin=1
      ipop=1
C
      icharge=0
C
      write1=.false.
      NORM=.true.
      SVD=.false.
      title='               '
c      natom=1
      nopt=0
      Nang=50
      DIRECT=.true.
      ATRHO=.false.
      SHFT=.false.
      ANG=.false.
      SHI=1.D0
      IDAMP=0
      EXTR=.false.
      istart=0
      SVD=.false.
      Exx=.true.
      Coul=.false.
      Scf1=.true.
      Prop=.false.
      iforce=0
      ngeo=0
      FTOL=1.0D-04
      delta = 0.01D0
      inmod = 0
      imin=0
      ibrent=0
      GRAD=.true.
      BSSE=.false.
      scale=1.0D0
      scale2=1.0D0
      thermo=0
      TEMP=298.0D0
      sigma=1
      popf=.false.
c
      h=0.15D0
      nsteps=10
      Nscale=20000
c
      sol=.false.
      free=.false.
      primera=.true.
      parsearch=.false.
      watermod=0
      noconverge=0
      converge=0
c
c      do n=1,ntq
c       map(n)=n
c      enddo
c calls generator of table for incomplete gamma functions
c
       call GENERF
       call GENERFS
       call GRIDLIO
       npas=0
c-----------------------------------------------------------------------
c reads input file
      
      if (.not.int_basis) then
      inquire(file=basis,exist=exists)
      if (.not.exists) then
      write(*,*) 'ERROR CANNOT FIND INPUT FILE ON UNIT 1',basis
      stop
      else
c  name of output file
      ikk=1
      do while (basis(ikk:ikk).ne.' ')
       ikk=ikk+1
      enddo
c
      ikk=ikk-1
c      name2=basis(1:ikk)//'.out'
c
      open(unit=1,file=basis,iostat=ios)
c      open(unit=2,file=output)
      endif
      endif

      open(unit=18,file=fcoord)
      open(unit=85,file=fmulliken)
      open(unit=88,file=frestart)
c-------------------------------------------------------
      date='date'
      write(*,*) 'JOB STARTED NOW'
      call system(date)
      do i=1,natom
        done(i)=.false.
      enddo

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
      if (.not.int_basis) then
        read(1,100) whatis
      endif
c      write(2,100) whatis
c
c
      No=0
      Nd=0
      Ne=0
      NBAS=0
      M=0
      Md=0

      allocate(natomc(natom),nnps(natom),nnpp(natom),nnp(natom))
      allocate(nnpd(natom),nns(natom),nnd(natom),atmin(natom))
      allocate(jatc(natom,natom))
      
      do i=1,natom
        natomc(i)=0 
      enddo
c-------------------------------------------------------------------------
c  BASIS SETS ------------------------------------------------------------
c-------------------------------------------------------------------------
      if (.not.int_basis) then
      do 25 while (whatis.ne.'endbasis')
c
        NBAS=NBAS+1
c signals if a basis set was not used
        used=.false.
        atmint=100000.
        read(1,*) iatom,nraw,ncon
c        write(2,600) iatom,nraw,ncon
c
c reads contraction scheme. The value for p,d ,f should not be repeated
c 3 ,6 , 10 .....   times
c reads also angular momentum for each of the contractions
c 0 for s 1 for p, etc
c
        read(1,*) (ncf(i),i=1,ncon)
c        write(2,*) (ncf(i),i=1,ncon)
        read(1,*) (lt(i),i=1,ncon)
c        write(2,*) (lt(i),i=1,ncon)
c
c loop over all primitives, no repeating p, d
        do i=1,nraw
          read(1,*) at(i),ct(i)
c       write(2,700) at(i),ct(i)
c       write(*,*) atmint,at(i)
          if(at(i).lt.atmint) atmint=at(i)
        enddo
c
        do j=1,natom
          if(Iz(j).eq.iatom.and.(.not.done(j))) then
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

c =====>>>>>>  M stores # of contractions  <<<<<<===========
            index=0
            do k=1,ncon
c
              M=M+Num(lt(k))

c nshell gives the # of functions s, p, d  etc
              nshell(lt(k))=nshell(lt(k))+Num(lt(k))
c
              do l2=1,Num(lt(k))
                No=No+1
c
c normalization
c
                if (NORM) then
                  do l=1,ncf(k)
c 
                    index=index+1

                    if(lt(k).eq.0) then
c
                      xnorm=sqrt((2.D0*at(index)/pi)**3)
                      xnorm=sqrt(xnorm)   
                      c(No,l)=ct(index)*xnorm
                      a(No,l)=at(index)
                    elseif(lt(k).eq.1) then
c
                      xnorm=sqrt((2.D0*at(index)/pi)**3)*4.D0*at(index)
                      xnorm=sqrt(xnorm)
                      c(No,l)=ct(index)*xnorm
                      a(No,l)=at(index)
                    elseif (lt(k).eq.2) then
c
                  xnorm=sqrt((2.D0*at(index)/pi)**3)*(4.D0*at(index))**2
                      xnorm=sqrt(xnorm)
                      c(No,l)=ct(index)*xnorm
                      a(No,l)=at(index)
                    endif
                  enddo
                else
c no normalization case
                  do l=1,ncf(k)
                    index=index+1
                    c(No,l)=ct(index)
                    a(No,l)=at(index)
                  enddo
                endif

c repeat the index only for p,d and f
                if (l2.ne.Num(lt(k))) then
                  index=index-ncf(k)
                endif
c
                Nuc(No)=j
                ncont(No)=ncf(k)
                nlb(No)=lt(k)
              enddo
            enddo
          endif
        enddo
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
c        write(2,*) iatom,nraw,ncon
c
c reads contraction scheme. The value for p,d ,f should not be repeated
c 3 ,6 , 10 .....   times. Reads also angular type , 
c 0 for s , 1 for p etc
        read(1,*) (ncf(i),i=1,ncon)
c        write(2,*) (ncf(i),i=1,ncon)
        read(1,*) (lt(i),i=1,ncon)
c        write(2,*) (lt(i),i=1,ncon)
c     
c
c loop over all primitives, repeating p, d
        do i=1,nraw
          read(1,*) at(i),ct(i)
c          write(2,700) at(i),ct(i)
        enddo
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
            Md=Md+Num(lt(k))
c          write(*,*) md 
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
 71               xnorm=sqrt((2.D0*at(index)/pi)**3)
                  xnorm=sqrt(xnorm)
                  cd(Nd,l)=ct(index)*xnorm
                  ad(Nd,l)=at(index)
c                  ad(Nd,l)=2.D0*at(index)
                  goto 47
c
 81               xnorm=sqrt((2.D0*at(index)/pi)**3)*4.D0*at(index)
                  xnorm=sqrt(xnorm)
                  cd(Nd,l)=ct(index)*xnorm
c                  ad(Nd,l)=2.D0*at(index)
                  ad(Nd,l)=at(index)
                  goto 47
c
 91               xnorm=sqrt((2.D0*at(index)/pi)**3)*(4.D0*at(index))**2
                  xnorm=sqrt(xnorm)
                  cd(Nd,l)=ct(index)*xnorm
c                  ad(Nd,l)=2.D0*at(index)
                  ad(Nd,l)=at(index)
                  goto 47
c
 47             continue
              else
c
c no normalization case
c
                do l=1,ncf(k)
                  index=index+1
c
                   cd(Nd,l)=ct(index)
c                  ad(Nd,l)=2.D0*at(index)
                  ad(Nd,l)=at(index)
                enddo
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
 46         continue
          endif
 45     continue
c
        read(1,100) whatis
c        write(2,100) whatis
 25   enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc  READING FROM INTERNAL LIO BASIS  cccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      else

      call getenv("LIOHOME",liohome)
      if (liohome == "") then
        write(*,*) "LIOHOME is not set! Cannot use basis set keywords"
        write(*,*) "Either set LIOHOME to your lio installation location
     > or specify a basis set file"
        stop
      endif
      int_basis_file=trim(liohome)//"/dat/basis/"//basis_set
      fit_basis_file=trim(liohome)//"/dat/basis/fitting/"//fitting_set
      !write(*,*) "INTERNAL BASIS FILE: ",trim(int_basis_file)
      !write(*,*) "INTERNAL FIT BASIS FILE: ",trim(fit_basis_file)

      inquire(file=int_basis_file,exist=exists)
      if (.not.exists) then
        write(*,*) 'THE BASIS SET ',trim(basis_set),
     >               ' COULD NOT BE FOUND'
        stop
      endif
      inquire(file=fit_basis_file,exist=exists)
      if (.not.exists) then
      write(*,*) 'THE FITTING BASIS SET ',trim(fitting_set),
     >             ' COULD NOT BE FOUND'
      stop
      endif

      open(unit=1,file=int_basis_file,iostat=ios)

      ! skip over all commented and blank lines
      inp_line = ""
      do while (inp_line == "")
        read(1,*) inp_line
        read(inp_line,'(a1)') inp_char
        if (inp_char == "#" .or. inp_line == "") then
          inp_line = ""
        endif
      enddo

      read(inp_line,100) whatis

      ! NOTE: we assume the entire section for an element (from "gaussian"
      ! to the last function values) is not broken by blank or commented lines
      do 26 while (whatis.ne.'endbasis')
c
        NBAS=NBAS+1
c signals if a basis set was not used
        used=.false.
        atmint=100000.
        read(1,*) iatom,nraw,ncon
        if (nraw.gt.nng) then
          write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          write(*,*) "This basis set contains an element with more",
     > " total primitives than your current maximum:",nng
          write(*,*) "Set nng in lioamber/liomods/garcha_mod.f to",
     > " a higher number, at least",nraw,"and recompile"
          write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          stop
        endif
c        write(2,600) iatom,nraw,ncon
c
c reads contraction scheme. The value for p,d ,f should not be repeated
c 3 ,6 , 10 .....   times
c reads also angular momentum for each of the contractions
c 0 for s 1 for p, etc
c
        read(1,*) (ncf(i),i=1,ncon)
c        write(2,*) (ncf(i),i=1,ncon)
        read(1,*) (lt(i),i=1,ncon)
c        write(2,*) (lt(i),i=1,ncon)
c
c loop over all primitives, no repeating p, d
        do i=1,nraw
          read(1,*) at(i),ct(i)
c       write(2,700) at(i),ct(i)
c       write(*,*) atmint,at(i)
          if(at(i).lt.atmint) atmint=at(i)
        enddo
c
        do j=1,natom
          if(Iz(j).eq.iatom.and.(.not.done(j))) then
          do i=1,ncon
            if (ncf(i).gt.nl) then
              write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
              write(*,*) "This basis set contains a function with more",
     > " primives than your currently set maximum: ",nl
              write(*,*) "Change nl in lioamber/liomods/param.f and",
     > " MAX_CONTRACTIONS in g2g/common.h to a larger number (at least",
     > ncf(i), ") and recompile"
              write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
              stop
            endif
          enddo
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

c =====>>>>>>  M stores # of contractions  <<<<<<===========
            index=0
            do k=1,ncon
c
              M=M+Num(lt(k))

c nshell gives the # of functions s, p, d  etc
              nshell(lt(k))=nshell(lt(k))+Num(lt(k))
c
              do l2=1,Num(lt(k))
                No=No+1
c
c normalization
c
                if (NORM) then
                  do l=1,ncf(k)
c 
                    index=index+1

                    if(lt(k).eq.0) then
c
                      xnorm=sqrt((2.D0*at(index)/pi)**3)
                      xnorm=sqrt(xnorm)   
                      c(No,l)=ct(index)*xnorm
                      a(No,l)=at(index)
                    elseif(lt(k).eq.1) then
c
                      xnorm=sqrt((2.D0*at(index)/pi)**3)*4.D0*at(index)
                      xnorm=sqrt(xnorm)
                      c(No,l)=ct(index)*xnorm
                      a(No,l)=at(index)
                    elseif (lt(k).eq.2) then
c
                  xnorm=sqrt((2.D0*at(index)/pi)**3)*(4.D0*at(index))**2
                      xnorm=sqrt(xnorm)
                      c(No,l)=ct(index)*xnorm
                      a(No,l)=at(index)
                    else
                write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                      write(*,*) "The basis set ",trim(basis_set),
     > " contains f (or higher) functions, and lio does not currently", 
     > " support them.  Choose another basis set."
                write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                      stop
                    endif
                  enddo
                else
c no normalization case
                  do l=1,ncf(k)
                    index=index+1
                    c(No,l)=ct(index)
                    a(No,l)=at(index)
                  enddo
                endif

c repeat the index only for p,d and f
                if (l2.ne.Num(lt(k))) then
                  index=index-ncf(k)
                endif
c
                Nuc(No)=j
                ncont(No)=ncf(k)
                nlb(No)=lt(k)
              enddo
            enddo
          endif
        enddo
c     
        if (.not.used.and.VERBOSE) then
          write(*,200) iatom
        endif

        ! skip over all commented and blank lines
        inp_line = ""
        do while (inp_line == "")
          read(1,*) inp_line
          read(inp_line,'(a1)') inp_char
          if (inp_char == "#" .or. inp_line == "") then
            inp_line = ""
          endif
        enddo
        read(inp_line,100) whatis

26    enddo

      close(1)

      open(unit=1,file=fit_basis_file,iostat=ios)

      ! skip over all commented and blank lines
      inp_line = ""
      do while (inp_line == "")
        read(1,*) inp_line
        read(inp_line,'(a1)') inp_char
        if (inp_char == "#" .or. inp_line == "") then
          inp_line = ""
        endif
      enddo

      read(inp_line,100) whatis
c
c Exactly the same should be repeated for charge density
c and exchange correlation potential basis sets
c
c CHARGE DENSITY --------------------------------------------------
c
      ! NOTE: we don't check right now for total primitives in the cd basis
      ! being less than nng...probably not a problem...
      do 27 while (whatis.ne.'endbasis')
        read(1,*) iatom,nraw,ncon
c        write(2,*) iatom,nraw,ncon
c
c reads contraction scheme. The value for p,d ,f should not be repeated
c 3 ,6 , 10 .....   times. Reads also angular type , 
c 0 for s , 1 for p etc
        read(1,*) (ncf(i),i=1,ncon)
c        write(2,*) (ncf(i),i=1,ncon)
        read(1,*) (lt(i),i=1,ncon)
c        write(2,*) (lt(i),i=1,ncon)
c     
c
c loop over all primitives, repeating p, d
        do i=1,nraw
          read(1,*) at(i),ct(i)
c          write(2,700) at(i),ct(i)
        enddo
c
        do 48 j=1,natom
          if (Iz(j).eq.iatom) then
c
c Mdd stores # of contractions in final basis, counting all possibilities
c for p , d etc
c
          index=0
          do 49 k=1,ncon
c
            Md=Md+Num(lt(k))
c          write(*,*) md 
            nshelld(lt(k))=nshelld(lt(k))+Num(lt(k))
c
            do 49 l2=1,Num(lt(k))
c
              Nd=Nd+1
c
              if (NORM) then
                do 50 l=1,ncf(k)
                  index=index+1
c
                  goto (72,82,92) lt(k)+1
                write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                      write(*,*) "The basis set ",trim(fitting_set),
     > " contains f (or higher) functions, and lio does not currently", 
     > " support them.  Choose another basis set"
                write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                  stop
c
 72               xnorm=sqrt((2.D0*at(index)/pi)**3)
                  xnorm=sqrt(xnorm)
                  cd(Nd,l)=ct(index)*xnorm
                  ad(Nd,l)=at(index)
c                  ad(Nd,l)=2.D0*at(index)
                  goto 50
c
 82               xnorm=sqrt((2.D0*at(index)/pi)**3)*4.D0*at(index)
                  xnorm=sqrt(xnorm)
                  cd(Nd,l)=ct(index)*xnorm
c                  ad(Nd,l)=2.D0*at(index)
                  ad(Nd,l)=at(index)
                  goto 50
c
 92               xnorm=sqrt((2.D0*at(index)/pi)**3)*(4.D0*at(index))**2
                  xnorm=sqrt(xnorm)
                  cd(Nd,l)=ct(index)*xnorm
c                  ad(Nd,l)=2.D0*at(index)
                  ad(Nd,l)=at(index)
                  goto 50
c
 50             continue
              else
c
c no normalization case
c
                do l=1,ncf(k)
                  index=index+1
c
                   cd(Nd,l)=ct(index)
c                  ad(Nd,l)=2.D0*at(index)
                  ad(Nd,l)=at(index)
                enddo
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
 49         continue
          endif
 48     continue
c
        ! skip over all commented and blank lines
        inp_line = ""
        do while (inp_line == "")
          read(1,*) inp_line
          read(inp_line,'(a1)') inp_char
          if (inp_char == "#" .or. inp_line == "") then
            inp_line = ""
          endif
        enddo

        read(inp_line,100) whatis
c        write(2,100) whatis
 27   enddo
      close(1)
      endif
c----- DIMENSION CONTROLS ------------------------------------
c
      iprob=0
      if (M.gt.ngDyn.or.M.gt.ng) then
        write(*,*) 'DIMENSION PROBLEMS WITH BASIS SET SIZE PARAMETER NG'
        write(*,*) 'NUMBER BASIS FUNCTIONS =',M,'<',ngDyn
        iprob=1
      endif
c
      if (Md.gt.ngdDyn.or.Md.gt.ngd) then
       write(*,*) 'DIMENSION PROBLEMS WITH AUXILIAR BASIS PARAMETER NGD'
        write(*,*) 'NUMBER AUXILIAR BASIS FUNCTIONS =',Md,'<',ngdDyn,ngd
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
        l=nlb(i)+1
        goto (11,22,33) l
c
c s case
  11    continue
c   
        Nucx(is)=Nuc(i)
        indexii(is)=i
        ncontx(is)=ncont(i)
c  
        do j=1,ncontx(is)
          cx(is,j)=c(i,j)
          ax(is,j)=a(i,j)
        enddo
c
        is=is+1
        goto 44   
c
c p case
  22    continue
c   
        Nucx(ip)=Nuc(i)
        indexii(ip)=i
        ncontx(ip)=ncont(i)
c  
        do j=1,ncontx(ip)
          cx(ip,j)=c(i,j)
          ax(ip,j)=a(i,j)
        enddo
c
        ip=ip+1
        goto 44
c
c d case
  33    continue
c   
        Nucx(id)=Nuc(i)
        indexii(id)=i

        ncontx(id)=ncont(i)
c  
        do j=1,ncontx(id)
c
          cx(id,j)=c(i,j)
          ax(id,j)=a(i,j)
        enddo
c
        id=id+1
        goto 44
c
  44    continue
 210  continue
c
c final normalization for d
c 3 factor for x^2, y^2 and z^2
c is in the integral part, so shell structure could be used
c
c
c Now, it has to put temporary things into the real one
c
      do i=1,M
        Nuc(i)=Nucx(i)
        ncont(i)=ncontx(i)
c
        do j=1,ncont(i)
          c(i,j)=cx(i,j)
          a(i,j)=ax(i,j)
        enddo
      enddo

c same process, but for electronic density basis set ------------
c
      is=1
      ip=nshelld(0)+1
      id=nshelld(0)+nshelld(1)+1
c
c loop over the total # of basis
      do 310 i=1,Md
c
        l=nld(i)+1
        goto (111,222,333) l
c
c s case
 111    continue
c   
        Nucx(is)=Nucd(i)
        indexiid(is)=i
        ncontx(is)=ncontd(i)
c  
        do j=1,ncontx(is)
          cx(is,j)=cd(i,j)
          ax(is,j)=ad(i,j)
        enddo
c
        is=is+1
        goto 444   
c
c p case
 222    continue
c   
        Nucx(ip)=Nucd(i)
        indexiid(ip)=i
        ncontx(ip)=ncontd(i)
c  
        do j=1,ncontx(ip)
          cx(ip,j)=cd(i,j)
          ax(ip,j)=ad(i,j)
        enddo
c
        ip=ip+1
        goto 444
c
c d case
 333    continue
c   
        Nucx(id)=Nucd(i)
        indexiid(id)=i
        ncontx(id)=ncontd(i)
c  
        do j=1,ncontx(id)
          cx(id,j)=cd(i,j)
          ax(id,j)=ad(i,j)
        enddo
c
        id=id+1
        goto 444
c
 444    continue
 310  continue
c
c final normalization for d
c 3 factor for x^2, y^2 and z^2
c is in the integral part, so shell structure could be used
c
c
c Now, it has to put temporary things into the real one
c
      do i=1,Md
        Nucd(i)=Nucx(i)
        ncontd(i)=ncontx(i)
c
        do j=1,ncontd(i)
          cd(i,j)=cx(i,j)
          ad(i,j)=ax(i,j)
        enddo
      enddo
c
c---------------------------------------------------------
c POINTERS -----------------------------------------------
c
c first P
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2

c now F alpha
      M3=MM+1
c now S, F beta also uses the same position after S was used
      M5=MM+M3
c now G
      M7=MM+M5
c now Gm
      M9=M7+MMd
c now H
      M11=M9+MMd
c W ( eigenvalues ), also this space is used in least squares
      M13=MM+M11
c aux ( vector for ESSl)
      M15=M+M13
c Least squares
      M17=MM+M15
c vectors of MO alpha
      M18=MMd+M17
c vectors of MO beta
c     M18b=M18+M*NCOa
c------------------------------------------------------------
c
c Density matrix  construction - For closed shell only <<<<=========
c
c      if (Scf1) then
c      read(1,nml=scfinp)
c      endif
c
c variables defined in namelist cannot be in common ?
      OPEN1=OPEN
c      NCO2=NCO
      Nunp2=Nunp
c
      idip1=idip
      ipop1=ipop
      icharge1=icharge
      ispin1=ispin
c      do n=1,ntq
c        map1(n)=map(n)
c      enddo
c
c      if (ANG) then
c        a01=a0/.529177D0
c      else
c        a01=a0
c      endif
       
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
c
c DIMENSION TESTS -----------------------------------------------
c
      Ndim=5*M*(M+1)/2+3*Md*(Md+1)/2+M+M*NCO!+M*Ngrid
      if(verbose) write(*,*) 'en drive', M,Md,NCO
      if (Ndim.gt.ng2) then
        write(*,*) 'DIMENSION PROBLEMS WITH DYNAMICAL VECTOR NG2',Ndim,ng2
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
        open(unit=89,file=frestartin)
c
        if (.not.OPEN) then
c reads vectors of MO coefficients, basis is in the same order as 
c the given in input
          do l=1,M
            read(89,*) (XX(l,n),n=1,NCO)
          enddo
c
c puts vectors in dynamical allocation (to be used later)
c
            kk=0
          do k=1,NCO
            do i=1,M
             kk=kk+1
             RMM(M18+kk-1)=XX(indexii(i),k)
            enddo
          enddo

          do i=1,M
            do j=1,M
              do k=1,NCO
                X(i,j)=X(i,j)+2.0D0*XX(i,k)*XX(j,k)
              enddo
            enddo
          enddo
c
c open shell case
        else
c
          NCOa=NCO
          NCOb=NCO+Nunp
          M18b=M18+M*NCOa
c alpha
          do l=1,M
            read(89,*) (XX(l,n),n=1,NCOa)
          enddo
c
          kk=M18-1
          do k=1,NCOa
            do i=1,M
              kk=kk+1
              RMM(kk)=XX(indexii(i),k)
            enddo
          enddo
c
c Density Matrix
c
          do 331 i=1,M
            do 331 j=1,M
              X(i,j)=0.0D0
c
              do 139 k=1,NCOa
                X(i,j)=X(i,j)+XX(i,k)*XX(j,k)
 139          continue
 331      continue
c
c beta
          do l=1,M
            read(89,*) (XX(l,n),n=1,NCOb)
          enddo
c
          kk=M18b-1
          do k=1,NCOb
            do i=1,M
              kk=kk+1
              RMM(kk)=XX(indexii(i),k)
            enddo
          enddo
c 
c Density Matrix
c
          do i=1,M
            do j=1,M
              do k=1,NCOb
                X(i,j)=X(i,j)+XX(i,k)*XX(j,k)
              enddo
            enddo
          enddo
c
        endif
      endif
c
c density matrix kept temporarily in S
c END of case in which density matrix is explicitly given in input
c--------------------------------------------------------------------
c
c case for initial guess constructed from atomic densities -------
      if (ATRHO) then
        k1=0
        l1=0
        NN=0
        do i=1,NBAS
c
c         read(1,nml=RHOINP)
c
          do j=1,nnat(i)
c
            do k=1,NAO
              k1=k1+1
c
              do l=1,NGF
                l1=NN+l
                kl=(k-1)*NGF+l
                XX(l1,k1)=atcoef(kl)
                oc2(k1)=OCC(k)
              enddo
            enddo
c
            NN=NN+NGF
          enddo
        enddo
c------------------------------------------------------------
c S used as scratch array here
        do i=1,M
          do j=1,M
            X(i,j)=0.
            do l=1,k1
              X(i,j)=X(i,j)+oc2(l)*XX(i,l)*XX(j,l)
            enddo 
          enddo
        enddo
c
c approximate vectors construction
c
        if (.not.OPEN) then
          kk=M18-1
          do k=1,NCO
            do i=1,M
              kk=kk+1
              RMM(kk)=XX(i,k) 
            enddo
          enddo
c
        else
c
          NCOa=NCO
          NCOb=NCO+Nunp
          M18b=M18+M*NCOa
c
          kk=M18-1
          do k=1,NCOa
            do i=1,M
              kk=kk+1
              RMM(kk)=XX(indexii(i),k)
            enddo
          enddo
c
          kk=M18b-1
          do k=1,NCOb
            do i=1,M
              kk=kk+1
              RMM(kk)=XX(indexii(i),k)
            enddo
          enddo 
        endif
      endif
c
c density matrix stored temporarily in S, then it should be changed
c according to the shell ordering of the basis set and kept in P
c------ end of option atomic densities -------------------------------------
c
c changes to the shell order ( s , p, d....)
      k=0
      do j=1,M
        do i=j,M
          k=k+1
          RMM(k)=X(indexii(i),indexii(j))
        enddo
      enddo
c
      k=0
      do j=1,M
        do i=j,M
          k=k+1
          if (i.ne.j) then
            RMM(k)=RMM(k)*2.D0
          endif
        enddo
      enddo
c
c---- reads exchange fit data -------------
c
      TMP1=ATRHO
      TMP2=VCINP
c
c      ATRHO=.FALSE.
C      VCINP=.TRUE.
c----------------------------
c      if (Scf1) then
c      write(2,nml=scfinp)
c      endif
c
c      if (exx) then
c      write(2,nml=EXCH)
c      endif
c
c      if (coul) then
c      write(2,nml=COULx)
c      endif
c
c      ATRHO=TMP1
c      VCINP=TMP2
c--------------------------
c
c------- G2G Initialization ---------------------
c
      ntqpru=natom
      ngpru=ng0*natom
c      write(*,*) 'estoooo',ngpru, ngDyn, ng0, natom

      if(OPEN) then
        allocate(rhoalpha(M*(M+1)/2),rhobeta(M*(M+1)/2))
      else
        allocate(rhoalpha(1),rhobeta(1))
      endif

      call g2g_parameter_init(NORM,natom,natom,ngDyn,!ngdDyn,
     >                        rqm,Rm2,Iz,Nr,Nr2,Nuc,
     >                        M,ncont,nshell,c,a,
     >                        RMM,M18,M5,M3,rhoalpha,rhobeta,
     >                        NCO,OPEN,Nunp,nopt,Iexch,
     >                        e_, e_2, e3, wang, wang2, wang3)

      call aint_parameter_init(Md, ncontd, nshelld, cd, ad, Nucd,
     >                         af, RMM, M9, M11, STR, FAC, rmax)

c      write(*,*) '======>>>> SALIENDO DE DRIVE <<<<=========='

c      if (parsearch) then
c        call g2g_reload_atom_positions(igrid2)
c        call g2g_solve_groups(1, E, 0)
c        goto 101
c        En el goto de arriba deberia hacer deallocate antes de irse
c      endif

c nopt 0 static SCF calculation --------------------------------------
c      if (nopt.eq.0) then

c      if (OPEN) then
c Nunp : number of unpaired electrons
c
c      call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
c     >         Nucd,Mdd,ncontd,nshelld,cd,ad,P,scratch,vec,E,
c     >   nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
c
c      else
c       write(*,*) 'llamo a SCF',natom,ntatom
c      call SCF()
c       write(*,*) 'llamé a SCF',natom,ntatom
c
c
c      endif
c      endif
c
c -- BSSE option -------------------------------------------------
c      if (nopt.eq.4) then
c 
c      write(*,*)
c      write(*,*) ' SCF ENERGY CALCULATION WITH GHOST ATOMS (BSSE)'
c
c it is convenient to use as a starting guess , the vectors from
c a full system calculation. The terms corresponding to the ghost
c atoms is set to 0
c
c
c      k=0
c      do j=1,M
c      do i=j,M
c       k=k+1
c       if (ighost(Nuc(i))*ighost(Nuc(j)).eq.0) then
c       P(k)=0.0D0
c       endif
c      enddo
c      enddo
c
c      BSSE=.true.
c      NCO2=nco1
c      Nunp2=Nunp1

c

c      if (OPEN1) then
c       OPEN=.true.
c      call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
c     >         Nucd,Mdd,ncontd,nshelld,cd,ad,P,scratch,vec,E,
c     >   nopt,OPEN,NMAX,NCO2,ATRHO,VCINP,SHFT,Nunp2,GOLD,told,write)
c
c      else
c      call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
c     >         Nucd,Mdd,ncontd,nshelld,cd,ad,P,scratch,vec,E,
c     >   nopt,OPEN,NMAX,NCO2,ATRHO,VCINP,SHFT,Nunp2,GOLD,told,write)
c
c      endif
c      endif
c
c------------------------------------------------------
c Standard (SCF) molecular dynamics option
c      if (nopt.eq.3) then
c      if (istart.eq.0) then
c       do kk=1,ntatom
c        vx(kk)=0.D0
c        vy(kk)=0.D0
c        vz(kk)=0.D0
c       enddo
c      endif
c
c      call MD2(MEMO,NORM,natom,Iz,r,v,Nuc,M,ncont,nshell,c,a,
c     >     Nucd,Mdd,ncontd,nshelld,cd,ad,P,scratch,vec,E,name1,ikk,
c     >   nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
c
c      endif
c
c--- Geometry optimization option
c
c      if (nopt.eq.2) then
c       call geom(MEMO,FTOL,imin,ibrent,NORM,natom,Iz,r,Nuc,M,ncont,
c     >   nshell,c,a,Nucd,Mdd,ncontd,nshelld,cd,ad,P,scratch,vec,E,Pm,
c     >           delta,inmod,name1,ikk,thermo,TEMP,sigma,
c     >   nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
c
c       endif
c
c-----------------------------------------------------------------------
c ---- Output results --------
c vector of molecular orbitals is going to be written in SCF subroutine
c
c      date='date'
c      write(*,*) 'JOB FINISHED'
c      call system(date)
c
c---------------------------------------------------
c---------------------------------------------------
       deallocate(X,XX)
       allocate(X(M,4*M),XX(Md,Md))
       allocate(old1(MM))

       allocate(old2(MM))
       allocate(old3(MM))

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
