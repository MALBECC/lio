! Convet a .pdb file to an hybrid2.0 input file
! using a force field defined in amber.parm file
! original version M. Marti
! modified version for hybrid2.0 Nick 2017
	program PDB2HYB
        implicit none

!actualizando variables, Nick
        logical :: existfile
	character(20) :: input, inpfile, outfile, coordfile, core !pdb2hyb input name
	integer :: nqm !number of QM residues
	integer :: natom  !number of QM atoms
        integer :: nlink !number of linkatoms
	integer :: natpdb !number of atoms in input .pdb file
	integer :: nspec !number of QM species
        integer :: numres_p !number of residues in force field file
	integer, dimension(:), allocatable :: resqm !number of residue that will be QM
	character*4, dimension(:), allocatable :: tty !part el residue that will be QM
	character*4, dimension(:), allocatable :: rname_p !nombre del residue en el campode fuersas
	integer, dimension(:), allocatable ::  atxres !cantidad de atomos en el residuo
	character :: formato*5 !formato del archivo de coordenadas
	character*2, dimension(:,:), allocatable :: attype_p
        character*4, dimension(:,:), allocatable :: atname_p
        integer :: coord_format !tipo de formato para las coordenadas
	double precision, dimension(:,:), allocatable :: coord
        integer, allocatable, dimension(:) :: indexqm

	character :: letra !auxiliar
        character*2 :: attype_paux !auxiliar
        character*4 :: atname_paux !auxiliar
	character :: titulo*20 !auxiliar
	integer :: i, i2,j,k,l !auxiliar
	integer :: maxatxres

!variable sde lectura de coordenadas
	character :: atom*4 !for check end of PDB
        character, allocatable, dimension(:) :: rname*4,atname*4 !residue and atom name
        integer, allocatable, dimension(:) :: atnu,resnu !residue and arom number
        double precision, allocatable, dimension(:,:) :: r !position
        character :: ch*1
	character :: new*4,old*4 !auxiliar for format change
	logical :: endpaso !auxiliar
        character :: ch2A*2,ch2B*2,ch4*4 !auxiliares para tipos de atomos
        logical :: search
	character*4 :: cresnum, cresnum2
	character :: CPa1 !for input of CP
	logical :: CP !CP=true remove extra atoms from C and N terminal
	logical :: remove !choose if remove atom from system

        character*2 :: nameout !nombre atomo QM
	integer :: Z !Z atomo QM
	logical, dimension(118) :: QMspec
        character*2, allocatable, dimension(:) :: at
	integer*2, allocatable, dimension(:) :: atsp
	integer :: io
	integer :: warnings
	character, dimension(:), allocatable :: attype*2
        character*2, allocatable, dimension(:) :: spec
        integer*2, allocatable, dimension(:) :: specz

        CHARACTER (LEN=3), DIMENSION(118) :: vec
        vec=(/'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ',  &
        'Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ','Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ',  &
        'Zr ','Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ',  &
        'Nd ','Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ',  &
        'Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ','Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ',  &
        'Fm ','Md ','No ','Lr ','Rf ','Db ','Sg ','Bh ','Hs ','Mm ','Ds ','Uuu','Uub','Uut','Uuq','Uup','Uuh','Uus','Uuo'/)


	CP=.false.

	call get_command_argument(1,input)
	call get_command_argument(2,CPa1)
	call logo()        
	
	warnings=0

!verificacion del campo de fuersas
       INQUIRE(FILE="amber.parm", EXIST=existfile)
       IF ( .NOT. existfile) THEN !verifica la existencia del archivo de parametros
        WRITE(*,*) "no se encuentra el campo de fuerzas"
        WRITE(*,*) "verifique que el archivo amber.parm"
        WRITE(*,*) "se encuentre en el directorio" 
        STOP
       END IF

       if (len_trim(input) .eq. 0) then
! LEctura SCRIPT FILE*******************************************
        write(*,*) 'Escriba el nombre del script a leer'
        write(*,*) ' EL script debe tener: '
        write(*,*) 'nombre del PDB , formato (VMD, hyb, amber)'
        write(*,*) 'numero de residuos que seran QM'
        write(*,*) 'los residuos por numero, typo GLYL,ALAL, etc'
        read(*,*) inpfile       
       else
	inpfile=input
       end if

	if (len_trim(CPa1) .ne.0) then
	  if (CPa1.eq."T" .or. CPa1.eq."t") CP=.true.
	  if (CPa1.eq."F" .or. CPa1.eq."f") CP=.false.
	  write(*,*) "cpvale", CP
	end if


!verificacion del SCRIPT
       INQUIRE(FILE=inpfile, EXIST=existfile)
       IF ( .NOT. existfile) THEN 
        WRITE(*,*) "no se encuentra el input"
        WRITE(*,*) "verifique que el archivo ", inpfile
        WRITE(*,*) "se encuentre en el directorio"
        STOP
       END IF

       open(unit=10,file='info.out')
	if (CP) write(10,*) "CP active, removing terminal atoms"

! Reding input*******************************************
       open(unit=100,file=inpfile)
       read(100,*) coordfile,formato !archivo con corrdenadas y formato
       read(100,*) nqm !number of QM residues
       allocate(resqm(nqm),tty(nqm))

       do i=1,nqm
         read(100,*) resqm(i),tty(i)
       enddo
       read(100,*) titulo,nlink
       close(100)

       write(10,*) 'INPUT file: ',inpfile,formato
       write(10,*) 'QM residues: ',(resqm(i),tty(i),i=1,nqm)


       INQUIRE(FILE=coordfile, EXIST=existfile)
       IF ( .NOT. existfile) THEN
         WRITE(*,*) "no se encuentra el archivo de coordenadas"
         WRITE(*,*) "verifique que el archivo ", coordfile
         WRITE(*,*) "se encuentre en el directorio"
         STOP
       END IF


! READING FORCE FIELD
	maxatxres=0
        open(unit=3,file='amber.parm')
        search=.true.
        do while(search)
         read(3,'(a8)') titulo
         if(titulo.eq.'residues') then
          read(3,*) numres_p
	  allocate(rname_p(numres_p), atxres(numres_p))
          do i=1,numres_p
            read(3,*) rname_p(i),atxres(i)
            do j=1,atxres(i)
              read(3,*) atname_paux,attype_paux
            enddo
	    if (atxres(i) .gt. maxatxres) maxatxres=atxres(i)
          enddo
          search=.false.
         endif
        enddo
	close(3) !ver de usar un rewind

	allocate(atname_p(numres_p,maxatxres),attype_p(numres_p,maxatxres))

        open(unit=3,file='amber.parm')
        search=.true.
        do while(search)
         read(3,'(a8)') titulo
         if(titulo.eq.'residues') then
          read(3,*) numres_p
          do i=1,numres_p
            read(3,*) rname_p(i),atxres(i)
            do j=1,atxres(i)
              read(3,*) atname_p(i,j),attype_p(i,j)
            enddo
          enddo
          search=.false.
         endif
        enddo
	close(3)
        write(10,*) 'End of Force Field read'


!  Crea el nombre del .fdf a partir del nombre del .pdb
	j=0
        do i=1,20
          letra=coordfile(i:i)
          if( (letra.eq.' ' .or. letra.eq.'.') .and. j.eq.0) then
            j=i-1
          endif
        enddo
        outfile(1:j)=coordfile(1:j)
        core(1:20)=' '
        core(1:j)=coordfile(1:j)
        outfile(j+1:j+4)='.fdf'

! abre el archivo de salida
        open(unit=2,file=outfile)


! ******************************************************************
! Lectura de coordenadas
	open(unit=1,file=coordfile)

	coord_format=0
	call check_format(coord_format,formato)
! coord_format = 1 VMD case
! coord_format = 2 HYB case

        if(coord_format.eq.1) then! caso formato tipo VMD (unico funcional por el momento)
          write(10,*) 'Reading ', coordfile,' in VMD format'
	  j=1
	  endpaso=.true.
	  do while (endpaso)
            read(1,'(A4)',iostat=io) atom
	    if (io.eq.0) then
	      j=j+1
              if(atom.eq.'END ') endpaso=.false.
              if(atom.eq.'TER ') j=j-1
	    else
	      endpaso=.false.
	      exit
	    end if
          enddo

          natpdb=j-1
	  close(1)
 
          write(10,*) 'number of atoms: ',natpdb

	  allocate(rname(natpdb+1), atname(natpdb+1), atnu(natpdb), resnu(natpdb+1), r(natpdb,3))
	  resnu=0

	  open(unit=1,file=coordfile)
          j=1
          endpaso=.true.
          do while (endpaso)
            read(1,'(A4,I7,1x,A4,1x,A4,A,A4,4x,3f8.3)',iostat=io)  &
                 atom,atnu(j),atname(j),rname(j),ch, &
                 cresnum,r(j,1),r(j,2),r(j,3)

            if (io.eq.0) then
              if (j.eq.1) then
                resnu(j)=1
                cresnum2=cresnum
              else
                if (cresnum .eq. cresnum2) then
                  resnu(j)=resnu(j-1)
                else
                  resnu(j)=resnu(j-1)+1
                  cresnum2=cresnum
                end if
              end if

              if(atom.eq.'END ') endpaso=.false.
              if(atom.eq.'TER ') j=j-1


	      call removelist(atname(j),rname(j),remove)

	      if (CP .and. remove) then
		warnings=warnings+1
		write(10,*) "Warning: removing atom",atnu(j),atname(j), &
               " from residue ", rname(j), &
               "for keep this atom use ./pdb2hyb pdb2hyb.in f"
		natpdb=natpdb-1
		j=j-1
	      end if
	      j=j+1
            else
              endpaso=.false.
              exit
            end if

	    if (resnu(j) .gt. 9999) then
!por ahora hybrid soporta hasta 9999 residuos
!para que soporte mas hay q cambiar la numeracion ya que el formato limita la numeracion a solo 4 digitos
		write(*,*) "your system have more than hybrid max residues"
		write(*,*) "reduce system to <= 9999 residues"
		STOP
	    end if
          enddo
 	  close(1)


        elseif(coord_format.eq.2) then 
!not tested by nick, this case is part of original code
!this will need some modifications
          write(10,*) 'Reading ', coordfile,' in HYB format'
          j=1
          endpaso=.true.
          do while (endpaso)
            read(1,'(A4,I7,1x,A4,1x,A4,A,I4,4x,3f8.3)',end=51) &
               atom,atnu(j),atname(j),rname(j),ch, &
               resnu(j),r(j,1),r(j,2),r(j,3)
            j=j+1
            if(atom.eq.'END ') endpaso=.false.
          enddo
 51       natpdb=j-1
          write(10,*) 'number of atoms: ',natpdb
         
        else 
          write(10,*) formato,' is not suported, use VMD or HYB'
	  STOP 'format is not suported, use VMD or HYB'
        endif
       
        write(10,*) 'End of coordinates read'
      
!********************************************************************************
! Transform VMD format to hybrid format
        if(coord_format.eq.1) then
	  write(10,*) "transformando VMD format to HYB format"
          do i=1,natpdb
            old=atname(i)
            ch=old(1:1)
            if(ch.eq.'1'.or.ch.eq.'2'.or. &
            ch.eq.'3') then
              if(old(3:3).eq.' ') then
                new(1:2)=old(2:3)
                new(3:3)=old(1:1)
              else 
                new(1:3)=old(2:4)
                new(4:4)=old(1:1)
              endif
              atname(i)=new

            elseif(ch.eq.' ') then
              ch=old(2:2)	
              if(ch.eq.'1'.or.ch.eq.'2'.or. &
              ch.eq.'3') then
                new(1:2)=old(3:4)
                new(3:3)=old(2:2)
                new(4:4)=' '
	      else
                new(1:3)=old(2:4)
                new(4:4)=' '
	      endif
              atname(i)=new
	    endif
	  enddo
        endif
!********************************************************************************


!********************************************************************************
! Defining QM subsystem
	allocate(indexqm(natpdb))
	write(10,*) "Defining QM subsystem"

        i2=1
        do i=1,nqm !QM residues
          do j=1,natpdb !total number of atoms
            if(resnu(j).eq.resqm(i)) then

              if(tty(i).eq.'ALL '.or.tty(i).eq.'ALL') then !si el residuo esta marcado como ALL todo el residuo es QM          
                indexqm(i2)=j
                i2=i2+1
              else
                do k=1,numres_p
                  if(rname_p(k).eq.tty(i)) then
                    do l=1,atxres(k)
! busca si el atomo en el residuo queda o no como clasico 
! o sea si esta en el residuo XXXL
                      if(atname(j).eq.atname_p(k,l)) then
                        rname(j)=tty(i)
                        goto 60
                      else
                      endif
! si ninguno de los atomso coicnide con el  del residuo XXXL    
! el atomo sea QM guarda indexqm
                    enddo
                    indexqm(i2)=j
                    i2=i2+1
                  endif
                enddo
              endif
            endif
! loop sobre atomos del pdb
  60      continue  
          enddo
        enddo
        natom=i2-1


! escribe los atomos QM en info.out
        write(10,*) 'QM atoms:'
        do i=1,natom
          write(10,*) atname(indexqm(i)), &
          rname(indexqm(i)),resnu(indexqm(i))
        enddo


!********************************************************************************
!  determiancion de tipos de atomos QM
	allocate(attype(natom+nlink))

        do i=1,natom!atomos totales QM
          attype(i)='  '
          do j=1,numres_p !residuos en el campo de fuersas
            if(rname(indexqm(i)).eq.rname_p(j)) then !verifica mismo coincidencia en el nombre de residue
              do k=1,atxres(j) !cantidad de atomos en el residuo
                if(atname(indexqm(i)).eq.atname_p(j,k)) then !verifica mismo coincidencia en el nombre del atomo
                  attype(i)=attype_p(j,k) !asigna el tpo de atomo del campo de fuersas
                  goto 70
                endif
              enddo
            endif
          enddo

          if(attype(i).eq.'  ') then
            write(10,*) 'WARNING: QM resiue: ',i,atname(indexqm(i)),rname(indexqm(i))
            write(10,*) 'not found in parm trying second pass'
	    warnings=warnings+1

            do j=1,numres_p
              ch4=rname(indexqm(i))
              ch2A=ch4(1:2)
              ch4=rname_p(j)
              ch2B=ch4(1:2)
              if(ch2A.eq.ch2B) then
                do k=1,atxres(j)
                  if(atname(indexqm(i)).eq.atname_p(j,k)) then
                    attype(i)=attype_p(j,k)
	            write(10,*) 'found coincidence for resname: ',ch2A, &
                    'in: ',rname_p(j), 'atom name: ', attype(i) 
                    goto 70
                  endif
                enddo
              endif
            enddo
          endif
        
          if(attype(i).eq.'  ') then
            write(*,*) 'QM resiue: ',i,atname(indexqm(i)),rname(indexqm(i))
            write(*,*) 'not found in parm file assigning default types'
            write(*,*) 'CHEK THEM'
	    Stop
          endif      
 70     continue
        enddo


!Pasa el Atname de los atomos QM a nombre del atomo
	QMspec=.false.
	allocate(at(natom), atsp(natom), coord(natom+nlink,3))

        do i=1,natom
          ch4=atname(indexqm(i))
	  call compare_atom_name(ch4, nameout, Z)
	  at(i)=nameout
	  write(*,*) "asigno: ", ch4, "a ", nameout
	  QMspec(Z)=.true. !lista de tipos de aomos QM
          do j=1,3
            coord(i,j)=r(indexqm(i),j)
          enddo
        enddo
        

!verifico que especies hay en el subsistema QM
	nspec=0
	do i=1,118
	  if (QMspec(i) .eq. .true.) then
	    nspec=nspec+1
	  end if
	end do

	allocate(spec(nspec), specz(nspec))

	j=0
	do i=1, 118
	  if (QMspec(i) .eq. .true.) then
	    j=j+1
	    spec(j)=vec(i)
	    specz(j)=i
          end if
	end do

	do i=1,natom
	  do j=1,nspec
	    if (at(i).eq.spec(j)) then
	      atsp(i) = j
	    endif
	  enddo
	enddo

         write(*,*) 'FIN DETRMINACION DE LAS ESPECIES'

! asighnacion coords y tipos de los linkatomos
! se fija que especie es H
        do i=1,nspec
          if(spec(i).eq.'H') k=i
        enddo

        do i=natom+1,natom+nlink
          do j=1,3
            coord(i,j)=0.0
          enddo
          attype(i)='HC'
          indexqm(i)=natpdb+1
          atsp(i)=k
        enddo

        atname(natpdb+1)=''
        rname(natpdb+1)='LINK'


!c ***************************************************************
!Writing .fdf file

        write (2,'(a15,a16)') 'SystemName     ',core 
        write (2,'(a15,a16)') 'SystemLabel    ',core    
	write (2,'(a20)')    'CenterMolecule T    '
        write (2,'(a20)')    'CG.Nick_center F    '
	write (2,'(a20,i3)') 'NumberOfAtoms       ',nlink+natom
	write (2,'(a20,i2)') 'NumberOfSpecies     ',nspec
	write (2,'(a20)')    'NetCharge 0.0       ' 	
	write (2,'(a20)')    'FixSpin  T          '
	write (2,'(a20)')    'TotalSpin 0.0       '
 	write (2,'(a30)')    'MD.MaxForceTol    0.04 eV/Ang  '	
	write (2,'(a22,i6)') 'NumberOfSolventAtoms  ',natpdb-natom

	if (nlink.ne.0) then
	  write (2,'(a20,i3)') 'LinkAtoms           ',nlink
	endif

	write (2,'(a30)')    '%block ChemicalSpeciesLabel    '

	do i=1,nspec
	  write (2,10) i,specz(i),spec(i) 
	enddo

 	write (2,'(a30)')    '%endblock ChemicalSpeciesLabel '
	Write (2,*)''
	write (2,'(a30)')    'AtomicCoordinatesFormat  Ang   '
	write (2,'(a40)')    '%block AtomicCoordinatesAndAtomicSpecies '

        do i=1,natom+nlink
	  write (2,'(3F10.5,i10)') coord(i,1),coord(i,2),coord(i,3),atsp(i)
	enddo

	write (2,'(a43)')    '%endblock AtomicCoordinatesAndAtomicSpecies' 
	write (2,*)''
	write (2,'(a25)')    '%block SoluteAtomTypes   '

        do i=1,natom+nlink
          write(2,'(2x,a2,2x,a4,2x,a4)') &
       attype(i),atname(indexqm(i)),rname(indexqm(i)) 
        enddo

	write (2,'(a25)')    '%endblock SoluteAtomTypes'
	write (2,*)''


        write (2,*)''
!        write (2,'(a26)')'PAO.BasisType split       '
!	write (2,'(a26)')'PAO.BasisSize       DZP   '
!        write (2,'(a27)')'PAO.EnergyShift      25 meV'
!        write (2,'(a26)')'PAO.SplitNorm       0.15  '
!        write (2,*)''
!        write (2,'(a26)')'WriteMDXmol F             '
!        write (2,'(a26)')'WriteMullikenPop 0        '
!        write (2,*)''
!	write (2,'(a26)')'XC.functional       GGA   '
!        write (2,'(a26)')'XC.authors          PBE   '
!        write (2,'(a27)')'MeshCutoff         150.0 Ry'
!        write (2,*)''
!        write (2,'(a26)')'SpinPolarized        T    '
!        write (2,'(a26)')'MaxSCFIterations    80    '
!        write (2,'(a26)')'DM.MixingWeight      0.1  '
!        write (2,'(a26)')'DM.NumberPulay       3    '
        write (2,'(a26)')'DM.UseSaveDM   T          '
        write (2,'(a26)')'MD.USeSaveXV   T          '
	write (2,'(a26)')'MD.UseSaveCG   T          '  
        write (2,'(a26)')'LongOutput    F           '
        write (2,*)''
        write (2,*)'OptimizationScheme 1'
        write (2,'(a26)')'MD.TypeOfRun    CG        '
        write (2,'(a26)')'MD.NumCGsteps   500       '
        write (2,*)''
	write (2,'(a26)')'%block CutOffRadius       '
	write (2,'(a26)')'QM  0.2                   '
	write (2,'(a26)')'QMMM 10.0                 '
	write (2,'(a26)')'MM 8.0                    '
	write (2,'(a26)')'BLO 8.0                   '
	write (2,'(a26)')'MMBLO 8.0                 '
	write (2,'(a26)')'%endblock CutOffRadius    '
	write (2,*)''
	
	write (2,*)''
	write (2,'(a22)')'%block SolventInput   '
        atom='ATOM'

	ch='H'
	k=1
        do j=1,natpdb
          do i=1,natom
            if(indexqm(i).eq.j) goto 80
          enddo
          write(2,'(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)') &
               atom,k,atname(j),rname(j),ch, &
               resnu(j),r(j,1),r(j,2),r(j,3)
	  k=k+1
 80     continue
        enddo
	write (2,'(a22)')'%endblock SolventInput'
	close(2)

	call write_lio_in()


	if (warnings.eq.0) then
	  write(*,*) 'PDB2HYB successfully finished'
	else
	  write(*,*) 'PDB2HYB finished'
	  write(*,*) 'there are: ', warnings, ' warnings'
	  write(*,*) 'check info.out'
	end if


	deallocate(resqm, tty, rname_p, atxres, atname_p, attype_p)
        deallocate(rname, atname, atnu, resnu, r, indexqm, attype)
        deallocate(at, atsp, coord, spec, specz)
	close(10)

   10   format(7x,i1,2x,i2,2x,A2) 	
	contains

	subroutine removelist(atname,rname,remove)
	implicit none
	character, intent(in) :: rname*4,atname*4
	logical, intent(out) :: remove
	remove=.false.
	if (atname .eq. " OXT") remove=.true.
	if (atname .eq. " H2" .or. atname .eq." H3") then
	  if (rname .ne. "WAT" .and. rname .ne. "HOH" .and. rname .ne. "PRE") remove=.true.
	end if 
	end subroutine removelist


	subroutine compare_atom_name(ch2, nameout, Z)
	implicit none
	CHARACTER (LEN=3), DIMENSION(118) :: vecCAP, vec
	character(len=2), intent(in) :: ch2
	character(len=2) :: ch2A
        character(len=2), intent(out) :: nameout
	integer :: i
        integer, intent(out) :: Z
	logical :: found, finish

        vecCAP=(/'H  ','He ','LI ','BE ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','MG ','AL ','SI ','P  ','S  ','CL ','AR ','K  ',  &
        'Ca ','SC ','TI ','V  ','CR ','MN ','FE ','CO ','NI ','CU ','ZN ','GA ','GE ','AS ','SE ','BR ','KR ','RB ','SR ','Y  ',  &
        'ZR ','Nb ','MO ','TC ','RU ','RH ','PD ','AG ','Cd ','IN ','SN ','SB ','TE ','I  ','XE ','CS ','BA ','LA ','Ce ','PR ',  &
        'Nd ','PM ','SM ','EU ','GD ','TB ','DY ','HO ','ER ','TM ','YB ','LU ','HF ','TA ','W  ','RE ','OS ','IR ','PT ','AU ',  &
        'Hg ','TL ','PB ','BI ','PO ','AT ','RN ','FR ','RA ','AC ','TH ','PA ','U  ','NP ','PU ','AM ','Cm ','BK ','CF ','ES ',  &
        'FM ','MD ','NO ','LR ','RF ','DB ','Sg ','BH ','HS ','MT ','DS ','UUU','UUB','UUT','UUQ','UUP','UUH','UUS','UUO'/)

        vec=(/'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ',  &
        'Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ','Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ',  &
        'Zr ','Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ',  &
        'Nd ','Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ',  &
        'Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ','Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ',  &
        'Fm ','Md ','No ','Lr ','Rf ','Db ','Sg ','Bh ','Hs ','Mm ','Ds ','Uuu','Uub','Uut','Uuq','Uup','Uuh','Uus','Uuo'/)

	i=0
	Z=0
	found=.false.
	finish=.false.

	ch2A=ch2

	do while (.not. finish)
	  i=i+1
	  if (ch2A .eq. vecCAP(i) .or. ch2A .eq. vec(i)) then
	    nameout=vec(i)
	    finish=.true.
	    found=.true.
	    Z=i
	  end if
	  if (i.gt.117) finish=.true.
	end do

	i=0
	if (.not. found) finish=.false.
	ch2A=ch2(1:1)
        do while (.not. finish)
          i=i+1
          if (ch2A .eq. vecCAP(i) .or. ch2A .eq. vec(i)) then
            nameout=vec(i)
            finish=.true.
            found=.true.
            Z=i
          end if
          if (i.gt.117) finish=.true.
        end do

	if (finish .and. .not. found) then
	  write(*,*) ch2, "not found"
	  STOP 
	end if
	end subroutine compare_atom_name


	subroutine check_format(vmd_corrd_format,formato)
        integer, intent(inout) :: vmd_corrd_format
        character, intent(in) :: formato*5
	coord_format=0
        if(formato.eq.'vmd') vmd_corrd_format=1
        if(formato.eq.'vmD') vmd_corrd_format=1
        if(formato.eq.'vMd') vmd_corrd_format=1
        if(formato.eq.'vMD') vmd_corrd_format=1
        if(formato.eq.'Vmd') vmd_corrd_format=1
        if(formato.eq.'VmD') vmd_corrd_format=1
        if(formato.eq.'VMd') vmd_corrd_format=1
        if(formato.eq.'VMD') vmd_corrd_format=1
        if(formato.eq.'hyb') vmd_corrd_format=2
        if(formato.eq.'hyB') vmd_corrd_format=2
        if(formato.eq.'hYb') vmd_corrd_format=2
        if(formato.eq.'hYB') vmd_corrd_format=2
        if(formato.eq.'Hyb') vmd_corrd_format=2
        if(formato.eq.'HyB') vmd_corrd_format=2
        if(formato.eq.'HYb') vmd_corrd_format=2
        if(formato.eq.'HYB') vmd_corrd_format=2
	end subroutine check_format

	subroutine write_lio_in()
	implicit none
	!agregar esto a futuro
	end subroutine write_lio_in
	


	subroutine logo
	implicit none
        write(*,*)
        write(*,1200)
        write(*,1201)
        write(*,1202)
        write(*,1203)
        write(*,1204)
        write(*,1205)
        write(*,*)
 1200 FORMAT(4x,"██████╗ ██████╗ ██████╗     ██████╗     ██╗  ██╗██╗   ██╗██████╗ ")
 1201 FORMAT(4x,"██╔══██╗██╔══██╗██╔══██╗    ╚════██╗    ██║  ██║╚██╗ ██╔╝██╔══██╗")
 1202 FORMAT(4x,"██████╔╝██║  ██║██████╔╝     █████╔╝    ███████║ ╚████╔╝ ██████╔╝")
 1203 FORMAT(4x,"██╔═══╝ ██║  ██║██╔══██╗    ██╔═══╝     ██╔══██║  ╚██╔╝  ██╔══██╗")
 1204 FORMAT(4x,"██║     ██████╔╝██████╔╝    ███████╗    ██║  ██║   ██║   ██████╔╝")
 1205 FORMAT(4x,"╚═╝     ╚═════╝ ╚═════╝     ╚══════╝    ╚═╝  ╚═╝   ╚═╝   ╚═════╝ ")
!ascii art page: http://patorjk.com/software/taag/#p=display&f=Ghost&t
	end subroutine logo

	end program


