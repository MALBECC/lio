!------------------------------------------------------------------------------!
! READS AND DRIVES SUBROUTINE                                                  !
! reads and prepare basis sets and geometry to perform electronic structure
! calculations using density functional theory and gaussian basis sets for
! the expansion of MO, charge density and exchange correlation potential
!
! Started 13 March 1992, Dario Estrin
!
!------------------------------------------------------------------------------!

      SUBROUTINE drive(ng2,ngDyn,ngdDyn)
      USE garcha_mod, ONLY : a,c, basis, done, done_fit, natomc, nnps,         &
      nnpp, nnpd, nns, nnp, nnd, atmin, jatc, ncf, lt, at, ct, nnat, nshell,   &
      nuc, ncont, nlb, nshelld, cd, ad, Nucd, ncontd, nld, Nucx, indexii,      &
      ncontx, cx, ax, indexiid, X, RMM, rhoalpha,rhobeta, af, charge,          &
      basis_set, fitting_set, e_, e_2, e3, exists, NORM, fcoord,               &
      fmulliken, natom, frestart, M, FAC, Iexch, int_basis, max_func, &
      frestartin, Md, NCO, nng, npas, Nr, used, STR, omit_bas, Nr2,   &
      wang, wang2, wang3, VCINP, OPEN, whatis, Num, Iz, pi,             &
      Rm2, rqm, rmax, Nunp, nl, nt, ng, ngd, restart_freq,             &
      writexyz, number_restr, restr_pairs,restr_index,restr_k,restr_w,restr_r0,&
      mulliken, MO_coef_at, MO_coef_at_b, use_libxc, ex_functional_id, &
      ec_functional_id

      USE ECP_mod, ONLY : ecpmode, asignacion
      USE fileio , ONLY : read_coef_restart, read_rho_restart
      use td_data    , only: td_do_pop
      use fileio_data, only: verbose, rst_dens
      use ghost_atoms_subs, only: summon_ghosts

      IMPLICIT NONE
      LOGICAL :: basis_check
      CHARACTER*255 :: int_basis_file, fit_basis_file
      CHARACTER*255 :: liohome
      CHARACTER*255 :: inp_line
      CHARACTER :: inp_char
      CHARACTER*3 :: simb

      INTEGER, INTENT(IN) :: ng2, ngDyn, ngdDyn
      REAL*8 :: atmint, iprob
      REAL*8 :: xnorm
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: restart_coef, restart_coef_b, &
                                             restart_dens, restart_adens,  &
                                             restart_bdens
      INTEGER :: NCOa, NCOb, ncon, nraw
      INTEGER :: is,ip,id, index
      INTEGER :: igpu, ios, NBAS, iatom
      INTEGER :: M3, M5, M7, M9, M11, M13, M15, M17, M18, M18b, MM, MMd !punteros
      INTEGER :: Nd, Ndim, nopt
      INTEGER :: i,ikk,j,k,k1,kk,kl,kkk,l,l1,l2,n,NN, No !auxiliars
!------------------------------------------------------------------------------!
! parameters for 2 basis sets, normal, density
! ng maximum # of contracted functions , nl maximum # of primitives in
! a contraction
! c, cd ce , linear coefficients
! a , ad, ae exponents
! Nuc , indicates center , ncont # of primitives in the correspondent
! contraction, and nx ny nz exponents of x , y and z in the cartesian
! gaussian
! ------------------
!
! Dimensions for Overlap and Fock matrices
!
!
! Angular momenta : up to f functions ( could be easily extended if
! necessary)
!
! --- Defaults ---------------------------------------------------------
!
! NORM true , expansion in normalized gaussians, so normalization factor
! included in coefficients
      nopt=0
!
! calls generator of table for incomplete gamma functions
!
       call GENERF
       call GENERFS
       call GRIDLIO
       npas=0
!-----------------------------------------------------------------------

! reads input file

      if (.not.int_basis) then
      inquire(file=basis,exist=exists)
      if (.not.exists) then
      write(*,*) 'ERROR CANNOT FIND INPUT FILE ON UNIT 1',basis
      stop
      else
!  name of output file
      ikk=1
      do while (basis(ikk:ikk).ne.' ')
       ikk=ikk+1
      enddo
!
      ikk=ikk-1
!      name2=basis(1:ikk)//'.out'
!
      open(unit=1,file=basis,iostat=ios)
!      open(unit=2,file=output)
      endif
      endif

      if (writexyz) open(unit=18,file=fcoord)
      if ((mulliken).or.(td_do_pop.gt.0)) open(unit=85,file=fmulliken)
      if (restart_freq.gt.0) open(unit=88,file=frestart)

!-------------------------------------------------------
      do i=1,natom
        done(i)=.false.
        done_fit(i)=.false.
      enddo

!c -------------------------------------------------------------
!c for each kind of atom, given by Z
!c
!c Basis set for MO expansion
!c
!c reads whatis: Gaussian for the time being,
!c      later stdbas will be included
!c reads iatom, atomic number , nraw number of primitive gaussians
!c , counting 3 p 6 d .., and ncon number of contractions counted in the
!c same way
!c if coefficients correspond to normalized gaussians
!c normalization factor will  multiply coefficients, so in  the program
!c expressions for non-normalized gaussians are used.
!c an option is included for using non-normalized gaussians
!c !The contractions don't need to be normalized, it will be done
!c automatically
!c
!c
      if (.not.int_basis) then
        read(1,100) whatis
      endif

      No=0
      Nd=0
      NBAS=0
      M=0
      Md=0

      allocate(natomc(natom),nnps(natom),nnpp(natom),nnp(natom))
      allocate(nnpd(natom),nns(natom),nnd(natom),atmin(natom))
      allocate(jatc(natom,natom))

      do i=1,natom
        natomc(i)=0
      enddo
!c-------------------------------------------------------------------------
!c  BASIS SETS ------------------------------------------------------------
!c-------------------------------------------------------------------------
      if (.not.int_basis) then
      do 25 while (whatis.ne.'endbasis')
!c
        NBAS=NBAS+1
!c signals if a basis set was not used
        used=.false.
        atmint=100000.
        read(1,*) iatom,nraw,ncon
!c        write(2,600) iatom,nraw,ncon
!c
!c reads contraction scheme. The value for p,d ,f should not be repeated
!c 3 ,6 , 10 .....   times
!c reads also angular momentum for each of the contractions
!c 0 for s 1 for p, etc
!c
        read(1,*) (ncf(i),i=1,ncon)
!c        write(2,*) (ncf(i),i=1,ncon)
        read(1,*) (lt(i),i=1,ncon)
!c        write(2,*) (lt(i),i=1,ncon)
!c
!c loop over all primitives, no repeating p, d
        do i=1,nraw
          read(1,*) at(i),ct(i)
          if(at(i).lt.atmint) atmint=at(i)
        enddo
!c
        do j=1,natom
          if(Iz(j).eq.iatom.and.(.not.done(j))) then
            nnat(NBAS)=nnat(NBAS)+1
            done(j)=.true.
            used=.true.
            atmin(j)=atmint

!c  cosas que puso nano para "atomos cerca"

            nns(j)=0
            nnp(j)=0
            nnd(j)=0

            do kkk=1,ncon
              if (lt(kkk).eq.0) nns(j)=nns(j)+Num(lt(kkk))
              if (lt(kkk).eq.1) nnp(j)=nnp(j)+Num(lt(kkk))
              if (lt(kkk).eq.2) nnd(j)=nnd(j)+Num(lt(kkk))
            enddo

!c      write(*,*) 'nns y etc',nns(j),nnp(j),nnd(j),nnps(j)
!c     > ,nnpp(j),nnpd(j)

!c =====>>>>>>  M stores # of contractions  <<<<<<===========
            index=0
            do k=1,ncon
!c
              M=M+Num(lt(k))

!c nshell gives the # of functions s, p, d  etc
              nshell(lt(k))=nshell(lt(k))+Num(lt(k))
!c
              do l2=1,Num(lt(k))
                No=No+1
!c
!c normalization
!c
                if (NORM) then
                  do l=1,ncf(k)
!c
                    index=index+1

                    if(lt(k).eq.0) then
!c
                      xnorm=sqrt((2.D0*at(index)/pi)**3)
                      xnorm=sqrt(xnorm)
                      c(No,l)=ct(index)*xnorm
                      a(No,l)=at(index)
                    elseif(lt(k).eq.1) then
!c
                      xnorm=sqrt((2.D0*at(index)/pi)**3)*4.D0*at(index)
                      xnorm=sqrt(xnorm)
                      c(No,l)=ct(index)*xnorm
                      a(No,l)=at(index)
                    elseif (lt(k).eq.2) then
!c
                  xnorm=sqrt((2.D0*at(index)/pi)**3)*(4.D0*at(index))**2
                      xnorm=sqrt(xnorm)
                      c(No,l)=ct(index)*xnorm
                      a(No,l)=at(index)
                    endif
                  enddo
                else
!c no normalization case
                  do l=1,ncf(k)
                    index=index+1
                    c(No,l)=ct(index)
                    a(No,l)=at(index)
                  enddo
                endif

!c repeat the index only for p,d and f
                if (l2.ne.Num(lt(k))) then
                  index=index-ncf(k)
                endif
!c
                Nuc(No)=j
                ncont(No)=ncf(k)
                nlb(No)=lt(k)
              enddo
            enddo
          endif
        enddo
!c
        if ( (.not.used) .and. (verbose.gt.4) ) then
          write(*,200) iatom
        endif
!c
!c Exactly the same should be repeated for charge density
!c and exchange correlation potential basis sets
!c
!c CHARGE DENSITY --------------------------------------------------
!c
        read(1,*) iatom,nraw,ncon
!c
!c reads contraction scheme. The value for p,d ,f should not be repeated
!c 3 ,6 , 10 .....   times. Reads also angular type ,
!c 0 for s , 1 for p etc
        read(1,*) (ncf(i),i=1,ncon)
        read(1,*) (lt(i),i=1,ncon)
!c
!c loop over all primitives, repeating p, d
        do i=1,nraw
          read(1,*) at(i),ct(i)
        enddo
!c
        do 45 j=1,natom
          if (Iz(j).eq.iatom) then
            done_fit(j)=.true.
!c
!c Mdd stores # of contractions in final basis, counting all possibilities
!c for p , d etc
!c
          index=0
          do 46 k=1,ncon
!c
            Md=Md+Num(lt(k))
            nshelld(lt(k))=nshelld(lt(k))+Num(lt(k))
!c
            do 46 l2=1,Num(lt(k))
!c
              Nd=Nd+1
!c
              if (NORM) then
                do 47 l=1,ncf(k)
                  index=index+1
!c
                  goto (71,81,91) lt(k)+1
!c
 71               xnorm=sqrt((2.D0*at(index)/pi)**3)
                  xnorm=sqrt(xnorm)
                  cd(Nd,l)=ct(index)*xnorm
                  ad(Nd,l)=at(index)
!c                  ad(Nd,l)=2.D0*at(index)
                  goto 47
!c
 81               xnorm=sqrt((2.D0*at(index)/pi)**3)*4.D0*at(index)
                  xnorm=sqrt(xnorm)
                  cd(Nd,l)=ct(index)*xnorm
!c                  ad(Nd,l)=2.D0*at(index)
                  ad(Nd,l)=at(index)
                  goto 47
!c
 91               xnorm=sqrt((2.D0*at(index)/pi)**3)*(4.D0*at(index))**2
                  xnorm=sqrt(xnorm)
                  cd(Nd,l)=ct(index)*xnorm
!c                  ad(Nd,l)=2.D0*at(index)
                  ad(Nd,l)=at(index)
                  goto 47
!c
 47             continue
              else
!c
!c no normalization case
!c
                do l=1,ncf(k)
                  index=index+1
!c
                   cd(Nd,l)=ct(index)
!c                  ad(Nd,l)=2.D0*at(index)
                  ad(Nd,l)=at(index)
                enddo
              endif
!c
!c repeat the index only for p,d and f, criterium l2<max(l2)
              if (l2.ne.Num(lt(k))) then
                index=index-ncf(k)
              endif
!c
              Nucd(Nd)=j
              ncontd(Nd)=ncf(k)
              nld(Nd)=lt(k)
!c
 46         continue
          endif
 45     continue
!c
        read(1,100) whatis
 25   enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccc  READING FROM INTERNAL LIO BASIS  cccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      else

      call getenv("LIOHOME",liohome)
      if (liohome == "") then
        write(*,*) "LIOHOME is not set! Cannot use basis set keywords"
        write(*,*) "Either set LIOHOME to your lio installation", &
       " location or specify a basis set file"
        stop
      endif
      int_basis_file=trim(liohome)//"/dat/basis/"//basis_set
      fit_basis_file=trim(liohome)//"/dat/basis/fitting/"//fitting_set
      !write(*,*) "INTERNAL BASIS FILE: ",trim(int_basis_file)
      !write(*,*) "INTERNAL FIT BASIS FILE: ",trim(fit_basis_file)

      inquire(file=int_basis_file,exist=exists)
      if (.not.exists) then
        write(*,*) 'THE BASIS SET ',trim(basis_set), &
      ' COULD NOT BE FOUND'
        stop
      endif
      inquire(file=fit_basis_file,exist=exists)
      if (.not.exists) then
      write(*,*) 'THE FITTING BASIS SET ',trim(fitting_set), &
      ' COULD NOT BE FOUND'
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
!c
        NBAS=NBAS+1
!c signals if a basis set was not used
        used=.false.
        atmint=100000.
        read(1,*) iatom,nraw,ncon
        if (nraw.gt.max_func) then         !aca hay algo raro, Nick usa nraw!!!!
          write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          write(*,*) "This basis set contains an element with more", &
       " total primitives than your current maximum:",nng

	write(*,*) "iatom vale ", iatom

          write(*,*) "Set nng in lioamber/liomods/garcha_mod.f to", &
       " a higher number, at least",nraw,"and recompile"
          write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          stop
        endif
!c
!c reads contraction scheme. The value for p,d ,f should not be repeated
!c 3 ,6 , 10 .....   times
!c reads also angular momentum for each of the contractions
!c 0 for s 1 for p, etc
!c
        read(1,*) (ncf(i),i=1,ncon)
        read(1,*) (lt(i),i=1,ncon)
!c
!c loop over all primitives, no repeating p, d
        do i=1,nraw
          read(1,*) at(i),ct(i)
          if(at(i).lt.atmint) atmint=at(i)
        enddo
!c
        do j=1,natom
          if(Iz(j).eq.iatom.and.(.not.done(j))) then
          do i=1,ncon
            if (ncf(i).gt.nl) then
              write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
              write(*,*) "This basis set contains a function with ",&
      "more primives than your currently set maximum: ",nl
              write(*,*) "Change nl in lioamber/liomods/param.f and",&
      " MAX_CONTRACTIONS in g2g/common.h to a larger number (at least",&
      ncf(i), ") and recompile"
              write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
              stop
            endif
          enddo
            nnat(NBAS)=nnat(NBAS)+1
            done(j)=.true.
            used=.true.
            atmin(j)=atmint

!c  cosas que puso nano para "atomos cerca"

            nns(j)=0
            nnp(j)=0
            nnd(j)=0

            do kkk=1,ncon
              if (lt(kkk).eq.0) nns(j)=nns(j)+Num(lt(kkk))
              if (lt(kkk).eq.1) nnp(j)=nnp(j)+Num(lt(kkk))
              if (lt(kkk).eq.2) nnd(j)=nnd(j)+Num(lt(kkk))
            enddo

!c =====>>>>>>  M stores # of contractions  <<<<<<===========
            index=0
            do k=1,ncon
!c
              M=M+Num(lt(k))

!c nshell gives the # of functions s, p, d  etc
              nshell(lt(k))=nshell(lt(k))+Num(lt(k))
!c
              do l2=1,Num(lt(k))
                No=No+1
!c
!c normalization
!c
                if (NORM) then
                  do l=1,ncf(k)
!c
                    index=index+1

                    if(lt(k).eq.0) then
!c
                      xnorm=sqrt((2.D0*at(index)/pi)**3)
                      xnorm=sqrt(xnorm)
                      c(No,l)=ct(index)*xnorm
                      a(No,l)=at(index)
                    elseif(lt(k).eq.1) then
!c
                      xnorm=sqrt((2.D0*at(index)/pi)**3)*4.D0*at(index)
                      xnorm=sqrt(xnorm)
                      c(No,l)=ct(index)*xnorm
                      a(No,l)=at(index)
                    elseif (lt(k).eq.2) then
!c
                  xnorm=sqrt((2.D0*at(index)/pi)**3)*(4.D0*at(index))**2
                      xnorm=sqrt(xnorm)
                      c(No,l)=ct(index)*xnorm
                      a(No,l)=at(index)
                    else
                write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                      write(*,*) "The basis set ",trim(basis_set), &
      " contains f (or higher) functions, and lio does not currently", &
      " support them.  Choose another basis set."
                write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                      stop
                    endif
                  enddo
                else
!c no normalization case
                  do l=1,ncf(k)
                    index=index+1
                    c(No,l)=ct(index)
                    a(No,l)=at(index)
                  enddo
                endif

!c repeat the index only for p,d and f
                if (l2.ne.Num(lt(k))) then
                  index=index-ncf(k)
                endif
!c
                Nuc(No)=j
                ncont(No)=ncf(k)
                nlb(No)=lt(k)
              enddo
            enddo
          endif
        enddo
!c
        if (.not.used.and.(verbose.gt.4).and. .not.omit_bas) then
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
!c
!c Exactly the same should be repeated for charge density
!c and exchange correlation potential basis sets
!c
!c CHARGE DENSITY --------------------------------------------------
!c
      ! NOTE: we don't check right now for total primitives in the cd basis
      ! being less than nng...probably not a problem...
      do 27 while (whatis.ne.'endbasis')
        read(1,*) iatom,nraw,ncon
!c        write(2,*) iatom,nraw,ncon
!c
!c reads contraction scheme. The value for p,d ,f should not be repeated
!c 3 ,6 , 10 .....   times. Reads also angular type ,
!c 0 for s , 1 for p etc
        read(1,*) (ncf(i),i=1,ncon)
        read(1,*) (lt(i),i=1,ncon)
!c
!c loop over all primitives, repeating p, d
        do i=1,nraw
          read(1,*) at(i),ct(i)
        enddo
!c
        do 48 j=1,natom
          if (Iz(j).eq.iatom) then
             done_fit(j)=.true.
!c
!c Mdd stores # of contractions in final basis, counting all possibilities
!c for p , d etc
!c
          index=0
          do 49 k=1,ncon
!c
            Md=Md+Num(lt(k))
            nshelld(lt(k))=nshelld(lt(k))+Num(lt(k))
!c
            do 49 l2=1,Num(lt(k))
!c
              Nd=Nd+1
!c
              if (NORM) then
                do 50 l=1,ncf(k)
                  index=index+1
!c
                  goto (72,82,92) lt(k)+1
                write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                      write(*,*) "The basis set ",trim(fitting_set), &
      " contains f (or higher) functions, and lio does not currently", &
      " support them.  Choose another basis set"
                write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                  stop
!c
 72               xnorm=sqrt((2.D0*at(index)/pi)**3)
                  xnorm=sqrt(xnorm)
                  cd(Nd,l)=ct(index)*xnorm
                  ad(Nd,l)=at(index)
                  goto 50
!c
 82               xnorm=sqrt((2.D0*at(index)/pi)**3)*4.D0*at(index)
                  xnorm=sqrt(xnorm)
                  cd(Nd,l)=ct(index)*xnorm
!c                  ad(Nd,l)=2.D0*at(index)
                  ad(Nd,l)=at(index)
                  goto 50
!c
 92               xnorm=sqrt((2.D0*at(index)/pi)**3)*(4.D0*at(index))**2
                  xnorm=sqrt(xnorm)
                  cd(Nd,l)=ct(index)*xnorm
!c                  ad(Nd,l)=2.D0*at(index)
                  ad(Nd,l)=at(index)
                  goto 50
!c
 50             continue
              else
!c
!c no normalization case
!c
                do l=1,ncf(k)
                  index=index+1
!c
                   cd(Nd,l)=ct(index)
!c                  ad(Nd,l)=2.D0*at(index)
                  ad(Nd,l)=at(index)
                enddo
              endif
!c
!c repeat the index only for p,d and f, criterium l2<max(l2)
              if (l2.ne.Num(lt(k))) then
                index=index-ncf(k)
              endif
!c
              Nucd(Nd)=j
              ncontd(Nd)=ncf(k)
              nld(Nd)=lt(k)
!c
 49         continue
          endif
 48     continue
!c
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
 27   enddo
      close(1)
      endif




!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      basis_check=.true.
      do i=1,natom
!c chekeo de lectura de base para todos los atomos, Nick
        if (.not. done(i)) then
           call asignacion(Iz(i),simb)
           write(*,*) simb," havent got a basis set, check basis file"
           basis_check=.false.
        end if
      enddo
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,natom
!c chekeo de lectura de base auxiliar para todos los atomos, Nick
        if (.not. done_fit(i)) then
           call asignacion(Iz(i),simb)
           write(*,*) simb," havent got a fitting basis set, check", &
      " auxiliar basis file"
           basis_check=.false.
        end if
      enddo
      if (.not. basis_check) stop
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!c----- DIMENSION CONTROLS ------------------------------------
!c
      iprob=0
      if (M.gt.ngDyn.or.M.gt.ng) then
        write(*,*) 'DIMENSION PROBLEMS WITH BASIS SET SIZE PARAMETER NG'
        write(*,*) 'NUMBER BASIS FUNCTIONS =',M,'<',ngDyn
        iprob=1
      endif
!c
      if (Md.gt.ngdDyn.or.Md.gt.ngd) then
       write(*,*) 'DIMENSION PROBLEMS WITH AUXILIAR BASIS PARAMETER NGD'
        write(*,*) 'NUMBER AUXILIAR BASIS FUNCTIONS =',Md,'<',ngdDyn,ngd
        iprob=1
      endif
!c
      if (natom.gt.nt) then
       write(*,*) 'DIMENSION PROBLEMS WITH NUMBER OF ATOMS PARAMETER NT'
        write(*,*) 'NUMBER OF ATOMS =',natom,'<',nt
        iprob=1
      endif
!c
!c -------------------------------------------------------------
!c -- Initial guess -------------------------------------------
!c
!c--- changes basis order to one in which the order is given by
!c the angular type :all s first, then all p, then all d ,......
!c also converts into normalized gaussians, including normalization
!c factor into linear coefficients
!c for d shell, in order to keep the same contraction for all the shell
!c x^2, y^2 and z^2 are normalized to 3, later in the integral part,
!c this should be considered
!c
      is=1
      ip=nshell(0)+1
      id=nshell(0)+nshell(1)+1
!c
!c standard basis set  ---------------------------------------------
!c
!c loop over the total # of basis
      do 210 i=1,M
        l=nlb(i)+1
        goto (11,22,33) l
!c
!c s case
  11    continue
!c
        Nucx(is)=Nuc(i)
        indexii(is)=i
        ncontx(is)=ncont(i)
!c
        do j=1,ncontx(is)
          cx(is,j)=c(i,j)
          ax(is,j)=a(i,j)
        enddo
!c
        is=is+1
        goto 44
!c
!c p case
  22    continue
!c
        Nucx(ip)=Nuc(i)
        indexii(ip)=i
        ncontx(ip)=ncont(i)
!c
        do j=1,ncontx(ip)
          cx(ip,j)=c(i,j)
          ax(ip,j)=a(i,j)
        enddo
!c
        ip=ip+1
        goto 44
!c
!c d case
  33    continue
!c
        Nucx(id)=Nuc(i)
        indexii(id)=i

        ncontx(id)=ncont(i)
!c
        do j=1,ncontx(id)
!c
          cx(id,j)=c(i,j)
          ax(id,j)=a(i,j)
        enddo
!c
        id=id+1
        goto 44
!c
  44    continue
 210  continue
!c
!c final normalization for d
!c 3 factor for x^2, y^2 and z^2
!c is in the integral part, so shell structure could be used
!c
!c
!c Now, it has to put temporary things into the real one
!c
      do i=1,M
        Nuc(i)=Nucx(i)
        ncont(i)=ncontx(i)
!c
        do j=1,ncont(i)
          c(i,j)=cx(i,j)
          a(i,j)=ax(i,j)
        enddo
      enddo

!c same process, but for electronic density basis set ------------
!
      is=1
      ip=nshelld(0)+1
      id=nshelld(0)+nshelld(1)+1
!c
!c loop over the total # of basis
      do 310 i=1,Md
!c
        l=nld(i)+1
        goto (111,222,333) l
!c
!c s case
 111    continue
!c
        Nucx(is)=Nucd(i)
        indexiid(is)=i
        ncontx(is)=ncontd(i)
!c
        do j=1,ncontx(is)
          cx(is,j)=cd(i,j)
          ax(is,j)=ad(i,j)
        enddo
!c
        is=is+1
        goto 444
!c
!c p case
 222    continue
!c
        Nucx(ip)=Nucd(i)
        indexiid(ip)=i
        ncontx(ip)=ncontd(i)
!c
        do j=1,ncontx(ip)
          cx(ip,j)=cd(i,j)
          ax(ip,j)=ad(i,j)
        enddo
!c
        ip=ip+1
        goto 444
!c
!c d case
 333    continue
!c
        Nucx(id)=Nucd(i)
        indexiid(id)=i
        ncontx(id)=ncontd(i)
!c
        do j=1,ncontx(id)
          cx(id,j)=cd(i,j)
          ax(id,j)=ad(i,j)
        enddo
!c
        id=id+1
        goto 444
!c
 444    continue
 310  continue
!c
!c final normalization for d
!c 3 factor for x^2, y^2 and z^2
!c is in the integral part, so shell structure could be used
!c
!c
!c Now, it has to put temporary things into the real one
!c
      do i=1,Md
        Nucd(i)=Nucx(i)
        ncontd(i)=ncontx(i)
!c
        do j=1,ncontd(i)
          cd(i,j)=cx(i,j)
          ad(i,j)=ax(i,j)
        enddo
      enddo
!c
!c---------------------------------------------------------
!c POINTERS -----------------------------------------------
!c
!c first P
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2

!c now F alpha
      M3=MM+1
!c now S, F beta also uses the same position after S was used
      M5=MM+M3
!c now G
      M7=MM+M5
!c now Gm
      M9=M7+MMd
!c now H
      M11=M9+MMd
!c W ( eigenvalues ), also this space is used in least squares
      M13=MM+M11
!c aux ( vector for ESSl)
      M15=M+M13
!c Least squares
      M17=MM+M15
!c vectors of MO alpha
      M18=MMd+M17
!c vectors of MO beta
!c
!c Density matrix  construction - For closed shell only <<<<=========
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (ecpmode) then !agregadas por Nick para lectura de ECP
           call lecturaECP()   !lee parametros
           CALL allocate_ECP() !allocatea la matriz de Fock de p-potenciales y el vector con los terminos de 1 electron sin corregir
           CALL ReasignZ() !reasigna las cargas de los nucleos removiendo la carga del core
        end if
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c DIMENSION TESTS -----------------------------------------------
!c
      Ndim=5*M*(M+1)/2+3*Md*(Md+1)/2+M+M*NCO!+M*Ngrid
      if (Ndim.gt.ng2) then
        write(*,*) 'DIMENSION PROBLEMS WITH DYNAMICAL VECTOR NG2',Ndim,ng2
        iprob=1
      endif
!c
      if (iprob.eq.1) then
         write(*,*) 'PAUSE IS A DELETED FEATURE'
!        pause
      endif

      ! Gets the number of occupied orbitals in a closed shell system (or
      ! Spin Up in an open shell system).
      call get_nco(Iz, natom, nco, charge, NUNP, OPEN)

      ! Allocates and initialises rhoalpha and rhobeta
      if(OPEN) then
        allocate(rhoalpha(M*(M+1)/2),rhobeta(M*(M+1)/2))
      else
        allocate(rhoalpha(1),rhobeta(1))
      endif
      rhoalpha(:) = 0.0d0
      rhobeta(:)  = 0.0d0


      ! Reads coefficient restart and builds density matrix. The MO
      ! coefficients are read in the same order as basis sets.
      ! Then vectors are put in dynamical allocation (to be used later)
      if (VCINP) then
         call g2g_timer_start('restart_read')

         open(unit=89, file=frestartin)
         if (rst_dens .gt. 0) then
            allocate(restart_dens(M, M))
            if (.not. OPEN) then
               call read_rho_restart(restart_dens, M, 89)
               call sprepack('L', M, RMM(1), restart_dens)
            else
               allocate(restart_adens(M,M), restart_bdens(M,M))
               call read_rho_restart(restart_adens, restart_bdens, M, 89)
               restart_dens = restart_adens + restart_bdens
               call sprepack('L', M, RMM(1)  , restart_dens)
               call sprepack('L', M, rhoalpha, restart_adens)
               call sprepack('L', M, rhobeta , restart_bdens)
               deallocate(restart_adens, restart_bdens)
            endif
            deallocate(restart_dens)
         else
            allocate(restart_dens(M, M), restart_coef(M, NCO))
            if (.not.OPEN) then
               call read_coef_restart(restart_coef, restart_dens, M, NCO, 89)

               kk = 0
               do k=1, NCO
               do i=1, M
                  kk = kk + 1
                  MO_coef_at(kk) = restart_coef(indexii(i), k)
               enddo
               enddo
            else
               NCOa = NCO
               NCOb = NCO + Nunp
               allocate(restart_coef_b(M, NCOb), restart_adens(M,M), &
                        restart_bdens(M,M))

               call read_coef_restart(restart_coef, restart_coef_b, &
                                      restart_dens, restart_adens,  &
                                      restart_bdens, M, NCOa, NCOb, 89)
               kk = 0
               do k=1, NCOa
               do i=1, M
                  kk = kk + 1
                  MO_coef_at(kk) = restart_coef(indexii(i), k)
               enddo
               enddo

               kk = 0
               do k=1, NCOb
               do i=1, M
                  kk = kk + 1
                  MO_coef_at_b(kk) = restart_coef_b(indexii(i), k)
               enddo
               enddo
               deallocate(restart_coef_b)
            endif

            ! Reorders by s, p, d.
            k = 0
            do j=1, M
            do i=j, M
               k = k + 1
               RMM(k)      = restart_dens(indexii(i), indexii(j))
               if (i.ne.j) then
                  RMM(k)      = RMM(k)*2.D0
               endif
            enddo
            enddo

            if (OPEN) then
               k = 0
               do j=1, M
               do i=j, M
                  k = k + 1
                  rhoalpha(k) = restart_adens(indexii(i), indexii(j))
                  rhobeta(k)  = restart_bdens(indexii(i), indexii(j))
                  if (i.ne.j) then
                     rhoalpha(k) = rhoalpha(k)*2.0D0
                     rhobeta(k)  = rhobeta(k)*2.0D0
                  endif
               enddo
               enddo
               deallocate(restart_adens, restart_bdens)
            endif
            deallocate(restart_dens, restart_coef)
         endif

         close(89)
         call g2g_timer_stop('restart_read')
      endif
      ! End of restart.

!c------- G2G Initialization ---------------------
!c
      call g2g_parameter_init(NORM,natom,natom,ngDyn, &
                             rqm,Rm2,Iz,Nr,Nr2,Nuc, &
                             M,ncont,nshell,c,a, &
                             RMM,M5,M3,rhoalpha,rhobeta, &
                             NCO,OPEN,Nunp,nopt,Iexch, &
                             e_, e_2, e3, wang, wang2, wang3, &
			                    use_libxc, ex_functional_id, ec_functional_id)
              call summon_ghosts(Iz, natom, verbose)

      call aint_query_gpu_level(igpu)
      if (igpu.gt.1) then
      call aint_parameter_init(Md, ncontd, nshelld, cd, ad, Nucd, &
      af, RMM, M9, M11, STR, FAC, rmax, Iz)
      endif
      allocate(X(M,4*M))


!--------------------------------------------------------------------------------------
      IF (number_restr.GT.0) THEN
! Distance Restrain parameters read
       ALLOCATE( restr_pairs(2,number_restr), restr_index(number_restr)&
       ,restr_k(number_restr), restr_w(number_restr),&
       restr_r0(number_restr))
        call read_restrain_params()
      END IF
!--------------------------------------------------------------------------------------




 100  format (A8)
 200  format ('  Basis set corresponding to Z = ', I3,' was not used.')
      return
      END SUBROUTINE drive




!This subroutine calculates the required dimension for many program variables
!ngnu max number of basis functions per atom
!ngdnu max number auxiliar of basis functions per atom
!max_func max number of total gausian functions per atom

      SUBROUTINE DIMdrive(ngnu,ngdnu)
       USE garcha_mod, ONLY: int_basis, natom, basis, whatis, exists, Iz, done &
       ,ncf, lt, Num, at, ct, done_fit, basis_set, fitting_set, nng, at, ct,   &
       max_func
       USE ECP_mod, ONLY : asignacion
       IMPLICIT NONE
       INTEGER, INTENT(OUT) :: ngnu,ngdnu
       INTEGER :: ios
       INTEGER :: iatom,nraw,ncon ! atomic number, number of primitive gaussians, number of contractions
       INTEGER :: i,j,k, M, Md !auxiliares

       CHARACTER*255 :: fit_basis_file, int_basis_file
       CHARACTER*255 :: liohome
       LOGICAL :: basis_check
       CHARACTER*255 :: inp_line
       CHARACTER :: inp_char
       CHARACTER*3 :: simb

       max_func=1
       ALLOCATE (ncf(1), lt(1))
       CALL reallocate_ncf_lt(1)

       DO i=1,natom
         done(i)=.false.
         done_fit(i)=.false.
       END DO



       IF (.NOT.int_basis) THEN   !case for provide a basis file
         INQUIRE (FILE=basis,EXIST=exists)
         IF (.NOT.exists) THEN
           WRITE(*,*) 'ERROR CANNOT FIND INPUT FILE ON UNIT 1',basis
           STOP
         ELSE
           OPEN(UNIT=1,FILE=basis,IOSTAT=ios)
         END IF
       END IF

!-------------------------------------------------------

! -------------------------------------------------------------
! Basis set for MO expansion
!
! reads whatis: Gaussian for the time being,
! reads iatom, atomic number , nraw number of primitive gaussians
! , counting 3 p 6 d .., and ncon number of contractions counted in the
! same way

      IF (.NOT.int_basis) THEN !case for provide a basis file
        READ(1,100) whatis
      END IF
      M=0
      Md=0

!-------------------------------------------------------------------------
!  BASIS SETS ------------------------------------------------------------
!-------------------------------------------------------------------------

      IF (.NOT.int_basis) THEN !case for provide a basis fil
        DO 25 WHILE (whatis.NE.'endbasis')
          READ(1,*) iatom,nraw,ncon  !

          IF (max_func .LT. nraw) max_func=nraw
          IF (ncon .GE. nng) call reallocate_ncf_lt(ncon)  !agregada para variables dinamicas

          READ(1,*) (ncf(i),i=1,ncon)! # of contractions of function i
          READ(1,*) (lt(i),i=1,ncon) !l of function i

          DO i=1,nraw
            READ(1,*)
          END DO

          DO j=1,natom
            IF(Iz(j) .EQ. iatom .AND.(.NOT.done(j))) THEN
              done(j)=.true.
              DO k=1,ncon
                M=M+Num(lt(k))
              END DO
            END IF
          END DO

! CHARGE DENSITY --------------------------------------------------
          READ(1,*) iatom,nraw,ncon

          IF (max_func .LT. nraw) max_func=nraw
          IF (ncon .GE. nng) call reallocate_ncf_lt(ncon) !agregada para variables dinamicas

          READ(1,*) (ncf(i),i=1,ncon)
          READ(1,*) (lt(i),i=1,ncon)
! loop over all primitives, repeating p, d
          DO i=1,nraw
            READ(1,*)
          END DO
!
          DO j=1,natom
            IF (Iz(j) .EQ. iatom) THEN
              done_fit(j)=.true.
              DO k=1,ncon
                Md=Md+Num(lt(k))
              END DO
            END IF
          END DO
          READ(1,100) whatis
 25     END DO
        CLOSE(1)

      ELSE

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccc  READING FROM INTERNAL LIO BASIS  cccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        CALL getenv("LIOHOME",liohome)
        IF (liohome == "") THEN
          WRITE(*,*) "LIOHOME is not set! Cannot use basis set keywords"
          WRITE(*,*) "Either set LIOHOME to your lio installation", &
          " location or specify a basis set file"
          STOP
        END IF

        int_basis_file=trim(liohome)//"/dat/basis/"//basis_set
        fit_basis_file=trim(liohome)//"/dat/basis/fitting/"//fitting_set
        !write(*,*) "INTERNAL BASIS FILE: ",trim(int_basis_file)
        !write(*,*) "INTERNAL FIT BASIS FILE: ",trim(fit_basis_file)

        INQUIRE(FILE=int_basis_file,EXIST=exists)
        IF (.NOT.exists) THEN
          WRITE(*,*) 'THE BASIS SET ',trim(basis_set), &
          ' COULD NOT BE FOUND'
          STOP
        ENDIF

        INQUIRE(FILE=fit_basis_file,EXIST=exists)
        IF (.NOT.exists) THEN
          WRITE(*,*) 'THE FITTING BASIS SET ',trim(fitting_set), &
          ' COULD NOT BE FOUND'
          STOP
        ENDIF
!
        OPEN(UNIT=1,FILE=int_basis_file,IOSTAT=ios)
!
! skip over all commented and blank lines
        inp_line = ""
        DO WHILE (inp_line == "")
          READ(1,*) inp_line
          READ(inp_line,'(a1)') inp_char
          IF (inp_char == "#" .OR. inp_line == "") THEN
            inp_line = ""
          END IF
        END DO
!
        READ(inp_line,100) whatis
!
        DO 26 WHILE (whatis.NE.'endbasis')
          READ(1,*) iatom,nraw,ncon

          IF (max_func .LT. nraw) max_func=nraw
          IF (ncon .GE. nng) call reallocate_ncf_lt(ncon) !agregada para variables dinamicas

          READ(1,*) (ncf(i),i=1,ncon)
          READ(1,*) (lt(i),i=1,ncon)

          DO i=1,nraw
            READ(1,*) !at(i),ct(i)
          END DO

          DO j=1,natom
            IF(Iz(j) .EQ. iatom .AND. (.NOT.done(j))) THEN
              done(j)=.true.
              DO k=1,ncon
                M=M+Num(lt(k)) ! =====>>>>>>  M stores # of contractions  <<<<<<===========
              END DO
            END IF
          END DO

          inp_line = ""
          DO WHILE (inp_line == "") ! skip over all commented and blank lines
            READ(1,*) inp_line
            READ(inp_line,'(a1)') inp_char
            IF (inp_char == "#" .OR. inp_line == "") THEN
              inp_line = ""
            END IF
          END DO
          READ(inp_line,100) whatis
26      END DO
        CLOSE(1)

!c
!c Exactly the same should be repeated for charge density
!c
!c CHARGE DENSITY --------------------------------------------------

        OPEN(UNIT=1,FILE=fit_basis_file,IOSTAT=ios)

        inp_line = ""
        DO WHILE (inp_line == "") ! skip over all commented and blank lines
          READ(1,*) inp_line
          READ(inp_line,'(a1)') inp_char
          IF (inp_char == "#" .OR. inp_line == "") THEN
            inp_line = ""
          END IF
        END DO
!
        READ(inp_line,100) whatis

        DO 27 WHILE (whatis .NE. 'endbasis')
          READ(1,*) iatom,nraw,ncon

          IF (max_func .LT. nraw) max_func=nraw
          IF (ncon .GE. nng) call reallocate_ncf_lt(ncon) !agregada para variables dinamicas

          READ(1,*) (ncf(i),i=1,ncon)
          READ(1,*) (lt(i),i=1,ncon)
          DO i=1,nraw
            READ(1,*) !at(i),ct(i)
          END DO
!c
          DO 48 j=1,natom
            IF (Iz(j).EQ.iatom) THEN
              done_fit(j)=.true.
              DO k=1,ncon
                Md=Md+Num(lt(k))
              END DO
	    END IF
 48       END DO

          inp_line = ""
          DO WHILE (inp_line == "") ! skip over all commented and blank lines
            READ(1,*) inp_line
            READ(inp_line,'(a1)') inp_char
            IF (inp_char == "#" .or. inp_line == "") THEN
              inp_line = ""
            END IF
          END DO
!
          READ(inp_line,100) whatis
 27     END DO
        CLOSE(1)
      END IF


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      basis_check=.true.
      DO i=1,natom ! chekeo de lectura de base para todos los atomos, Nick
        IF (.NOT. done(i)) THEN
           CALL asignacion(Iz(i),simb)
           WRITE(*,*) simb," havent got a basis set, check basis file"
           basis_check=.false.
        END IF
      END DO
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO i=1,natom ! chekeo de lectura de base auxiliar para todos los atomos, Nick
        IF (.NOT. done_fit(i)) THEN
           CALL asignacion(Iz(i),simb)
           WRITE(*,*) simb," havent got a fitting basis set, check", &
      " auxiliar basis file"
           basis_check=.false.
        END IF
      END DO
      IF (.NOT. basis_check) STOP
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ngnu=M
      ngdnu=Md
      ALLOCATE (at(max_func), ct(max_func))

 100  FORMAT (A8)
      RETURN
      END SUBROUTINE DIMdrive


      SUBROUTINE reallocate_ncf_lt(Dime)
       USE garcha_mod, ONLY: ncf, lt, nng
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: Dime
       nng=Dime+1
       DEALLOCATE (ncf, lt)
       ALLOCATE (ncf(nng), lt(nng))
      ENDSUBROUTINE reallocate_ncf_lt

subroutine get_nco(atom_Z, n_atoms, n_orbitals, n_unpaired, charge, open_shell)
   use ghost_atoms_subs, only: adjust_ghost_charge
   implicit none
   integer, intent(in)  :: n_atoms, n_unpaired, charge, atom_Z(n_atoms)
   logical, intent(in)  :: open_shell
   integer, intent(out) :: n_orbitals

   integer :: icount, nuc_charge, electrons

   nuc_charge = 0
   do icount = 1, n_atoms
      nuc_charge = nuc_charge + atom_Z(icount)
   enddo
   call adjust_ghost_charge(atom_Z, n_atoms, nuc_charge)

   electrons = nuc_charge - charge
   if ((.not.open_shell) .and. (mod(electrons,2).ne.0)) then
      write(*,'(A)') "  ERROR - DRIVE: Odd number of electrons in a "&
                     &"closed-shell calculation."
      write(*,'(A)') "  Please check system charge."
      stop
   endif

   n_orbitals = ((nuc_charge - charge) - n_unpaired)/2

   return
end subroutine get_nco
