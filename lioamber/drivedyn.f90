      SUBROUTINE dymdrive(ngnu,ngdnu)
       USE garcha_mod, ONLY: int_basis, natom, basis, whatis, exists, Iz, done, ncf, lt, Num, at, ct, done_fit, basis_set, fitting_set, nng, at, ct, max_func
       use ECP_mod, only : asignacion
       IMPLICIT NONE
       INTEGER :: ios
       INTEGER :: iatom,nraw,ncon ! atomic number, number of primitive gaussians, number of contractions
       INTEGER :: i,j,k !auxiliares
       INTEGER :: M, Md

       CHARACTER*255 :: fit_basis_file, int_basis_file
       CHARACTER*255 liohome

!      logical Exx, parsearch, basis_check
      LOGICAL :: basis_check
      character*255 inp_line
      character inp_char
      CHARACTER*3 simb
!

       INTEGER, INTENT(OUT) :: ngnu,ngdnu

      max_func=1
      ALLOCATE (ncf(1), lt(1))
      call reallocate_ncf_lt(1)

        write(*,*) "pase por aqui"

!c ------------------
!c parameters for 2 basis sets, normal, density 
!c ng maximum # of contracted functions , nl maximum # of primitives in
!c a contraction
!c c, cd ce , linear coefficients
!c a , ad, ae exponents
!c Nuc , indicates center , ncont # of primitives in the correspondent
!c contraction, and nx ny nz exponents of x , y and z in the cartesian
!c gaussian
!c ------------------
!c
!c Dimensions for Overlap and Fock matrices

!c Angular momenta : up to f functions ( could be easily extended if
!c necessary)
!c
!c-----------------------------------------------------------------------
!c reads input file
!      
      if (.not.int_basis) then
        inquire(file=basis,exist=exists)
        if (.not.exists) then
         write(*,*) 'ERROR CANNOT FIND INPUT FILE ON UNIT 1',basis
         stop
        else
         open(unit=1,file=basis,iostat=ios)
        endif
      endif
!
!c-------------------------------------------------------
      do i=1,natom
        done(i)=.false.
        done_fit(i)=.false.
      enddo
!
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
      if (.not.int_basis) then
        read(1,100) whatis
	write(*,*) "what1", whatis
      endif
      M=0
      Md=0
!
!c-------------------------------------------------------------------------
!c  BASIS SETS ------------------------------------------------------------
!c-------------------------------------------------------------------------
      if (.not.int_basis) then
        do 25 while (whatis.ne.'endbasis')
          read(1,*) iatom,nraw,ncon  !

          if (max_func .lt. nraw) max_func=nraw
          if (ncon .ge. nng) call reallocate_ncf_lt(ncon)  !agregada para variables dinamicas

          read(1,*) (ncf(i),i=1,ncon)! # of contractions of function i
          read(1,*) (lt(i),i=1,ncon) !l of function i
          do i=1,nraw
            read(1,*) 
          enddo

          do j=1,natom
            if(Iz(j).eq.iatom.and.(.not.done(j))) then
              done(j)=.true.
              do k=1,ncon
                M=M+Num(lt(k))
              enddo
            endif
          enddo

!c CHARGE DENSITY --------------------------------------------------
          read(1,*) iatom,nraw,ncon

          if (max_func .lt. nraw) max_func=nraw
          if (ncon .ge. nng) call reallocate_ncf_lt(ncon) !agregada para variables dinamicas

          read(1,*) (ncf(i),i=1,ncon)
          read(1,*) (lt(i),i=1,ncon)
!c loop over all primitives, repeating p, d
          do i=1,nraw
            read(1,*) 
          enddo
!c
          do  j=1,natom
            if (Iz(j).eq.iatom) then
              done_fit(j)=.true.
              do k=1,ncon
                Md=Md+Num(lt(k))
              end do
            endif
          end do
        read(1,100) whatis
 25   enddo
     close(1)


      else


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccc  READING FROM INTERNAL LIO BASIS  cccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	write(*,*) "entre a base interna"

      call getenv("LIOHOME",liohome)
      if (liohome == "") then
        write(*,*) "LIOHOME is not set! Cannot use basis set keywords"
        write(*,*) "Either set LIOHOME to your lio installation", & 
        " location or specify a basis set file"
        stop
      endif
      int_basis_file=trim(liohome)//"/dat/basis/"//basis_set
      fit_basis_file=trim(liohome)//"/dat/basis/fitting/"//fitting_set
!
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
!
      open(unit=1,file=int_basis_file,iostat=ios)
	write(*,*) "abri archivo de base"
!
! skip over all commented and blank lines
      inp_line = ""
      do while (inp_line == "")
        read(1,*) inp_line
        read(inp_line,'(a1)') inp_char
        if (inp_char == "#" .or. inp_line == "") then
          inp_line = ""
        endif
      enddo
!
      read(inp_line,100) whatis
!
      do 26 while (whatis.ne.'endbasis')
        read(1,*) iatom,nraw,ncon

        if (max_func .lt. nraw) max_func=nraw
        if (ncon .ge. nng) call reallocate_ncf_lt(ncon) !agregada para variables dinamicas
        write(*,*) iatom, ncon, nng
!	Write(*,*) "lei", iatom
!        if (nraw.gt.nng) then
!          write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!          write(*,*) "This basis set contains an element with more",
!     > " total primitives than your current maximum:",nng
!          write(*,*) "Set nng in lioamber/liomods/garcha_mod.f to",
!     > " a higher number, at least",nraw,"and recompile"
!          write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!          stop
!        endif
        read(1,*) (ncf(i),i=1,ncon)
        read(1,*) (lt(i),i=1,ncon)

!c
        do i=1,nraw
          read(1,*) !at(i),ct(i)
        enddo

        do j=1,natom
!		write(*,*) "llegue a j", j
!		write(*,*) Iz(j)   ! aca pincha por que no tiene Iz aun, arreglado!
          if(Iz(j).eq.iatom.and.(.not.done(j))) then
!	 	write(*,*) "pase el if"

              done(j)=.true.


!          do i=1,ncon
!            if (ncf(i).gt.nl) then
!              write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!              write(*,*) "This basis set contains a function with more",
!     > " primives than your currently set maximum: ",nl
!              write(*,*) "Change nl in lioamber/liomods/param.f and",
!     > " MAX_CONTRACTIONS in g2g/common.h to a larger number (at least",
!     > ncf(i), ") and recompile"
!              write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!              stop
!            endif
!          enddo
!            nnat(NBAS)=nnat(NBAS)+1
!
!c =====>>>>>>  M stores # of contractions  <<<<<<===========
            do k=1,ncon
              M=M+Num(lt(k))
            enddo
          endif
        enddo
!
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
!
!	write(*,*) "pase un do"
26    enddo
!
      close(1)


!c
!c Exactly the same should be repeated for charge density
!c and exchange correlation potential basis sets
!c
!c CHARGE DENSITY --------------------------------------------------

      open(unit=1,file=fit_basis_file,iostat=ios)
	write(*,*) "entre a base interna auxiliar"
! skip over all commented and blank lines
      inp_line = ""
      do while (inp_line == "")
        read(1,*) inp_line
        read(inp_line,'(a1)') inp_char
        if (inp_char == "#" .or. inp_line == "") then
          inp_line = ""
        endif
      enddo
!
      read(inp_line,100) whatis
!c
!      ! NOTE: we don't check right now for total primitives in the cd basis
!      ! being less than nng...probably not a problem...

      do 27 while (whatis.ne.'endbasis')
        read(1,*) iatom,nraw,ncon

        if (max_func .lt. nraw) max_func=nraw
        if (ncon .ge. nng) call reallocate_ncf_lt(ncon) !agregada para variables dinamicas

        read(1,*) (ncf(i),i=1,ncon)
        read(1,*) (lt(i),i=1,ncon)
        do i=1,nraw
          read(1,*) !at(i),ct(i)
        enddo
!c
        do 48 j=1,natom
          if (Iz(j).eq.iatom) then
             done_fit(j)=.true.
!c
          do 49 k=1,ncon
            Md=Md+Num(lt(k))
 49       end do
	  end if
 48     end do
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
!
        read(inp_line,100) whatis
!	write(*,*) "pase un do auxiliar"
 27   enddo
      close(1)
      endif
!
!
!
!
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
!
!
!
!c----- DIMENSION CONTROLS ------------------------------------
!c
!      iprob=0
!      if (M.gt.ngDyn.or.M.gt.ng) then
!        write(*,*) 'DIMENSION PROBLEMS WITH BASIS SET SIZE PARAMETER NG'
!        write(*,*) 'NUMBER BASIS FUNCTIONS =',M,'<',ngDyn
!        iprob=1
!      endif
!c
!      if (Md.gt.ngdDyn.or.Md.gt.ngd) then
!       write(*,*) 'DIMENSION PROBLEMS WITH AUXILIAR BASIS PARAMETER NGD'
!        write(*,*) 'NUMBER AUXILIAR BASIS FUNCTIONS =',Md,'<',ngdDyn,ngd
!        iprob=1
!      endif
!c
!      if (natom.gt.nt) then
!       write(*,*) 'DIMENSION PROBLEMS WITH NUMBER OF ATOMS PARAMETER NT'
!        write(*,*) 'NUMBER OF ATOMS =',natom,'<',nt
!        iprob=1
!      endif
!c
!c -------------------------------------------------------------
!c -- Initial guess -------------------------------------------
!
!c DIMENSION TESTS -----------------------------------------------
!c
!      Ndim=5*M*(M+1)/2+3*Md*(Md+1)/2+M+M*NCO!+M*Ngrid
!      if(verbose) write(*,*) 'en drive', M,Md,NCO
!      if (Ndim.gt.ng2) then
!        write(*,*) 'DIMENSION PROBLEMS WITH DYNAMICAL VECTOR NG2',Ndim,ng2
!        iprob=1
!      endif
!c
!      if (iprob.eq.1) then
!         write(*,*) 'PAUSE IS A DELETED FEATURE'
!         pause
!      endif
!
       ngnu=M
       ngdnu=Md
       allocate (at(max_func), ct(max_func)) 

        write(*,*) 
        write(*,*) 
        write(*,*) 
	write(*,*) "termine drivedyn", " M Md nng", M, Md, nng
        write(*,*)
        write(*,*)
        write(*,*)


 100  FORMAT (A8)
      RETURN
      END SUBROUTINE dymdrive


      SUBROUTINE reallocate_ncf_lt(Dime)
       USE garcha_mod, ONLY: ncf, lt, nng
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: Dime
       nng=Dime+1
       write(*,*) "cambie nng ", nng
       DEALLOCATE (ncf, lt)
       ALLOCATE (ncf(nng), lt(nng))
      ENDSUBROUTINE reallocate_ncf_lt
