!
! This module implements an interface to the FORTRAN logical unit
! system. Based on code by Richard Maine.
!
!
! Alberto Garcia, December 30, 1996
! Rewritten as a single subroutine 
! with multiple entry points, March 7, 1998
!
! This scheme is actually the closest in spirit to f90 modules, but
! in perfectly legal f77.
!
!---------------------------------------------------------------
!
! Actually put into a F90 module CJP 1/3/99
!
Module io
  !
  ! Logical unit management. Units 0 to min_lun-1 are "reserved",
  !  since most of the "typical" files (output, etc) use them.
  !
  ! Logical units min_lun to min_max are managed by this module.
  !
  Implicit None
  !
  ! Kind parameters 
  !
  Integer, Private, Parameter :: I4B = Selected_int_kind(9)
  Integer, Private, Parameter :: DP  = Kind(1.d0)
  Integer, Private, Parameter :: SP  = Kind(1.0)
  !
  !----------------------------------------------------------------
  !     Module variables
  !
  Integer(I4B), Private            :: stdout, stderr
  Integer(I4B), Private, Parameter :: min_lun=10,max_lun=99
  Integer(I4B), Private, Parameter :: nunits=max_lun-min_lun+1
  Logical, Private                 :: lun_is_free(min_lun:max_lun)
  !
  !-----------------------------------------------------------------
  !
  !     Internal and dummy variables
  !
  Integer(I4B), Private :: i, iostat
  Logical, Private      :: used, named, opened
  Character, Private    :: filename*50, form*11
  !
  !-----------------------------------------------------------------
  !     Initialization section
  !
  Data lun_is_free /nunits*.True./
  Data stdout, stderr /6,0/
  !-----------------------------------------------------------------
  !
  !     Executable routines
  !
  !     Simple interfaces to modify standard units
  !
Contains

  Subroutine io_seterr(unit)
    Integer(I4B), Intent(INOUT) :: unit
    stderr = unit
    Return
  End Subroutine io_seterr
  !
  Subroutine io_setout(unit)
    Integer(I4B), Intent(INOUT) :: unit
    stdout = unit
    Return
  End Subroutine io_setout
  !
  Subroutine io_geterr(unit)
    Integer(I4B), Intent(INOUT) :: unit
    unit = stderr
    Return
  End Subroutine io_geterr
  !
  Subroutine io_getout(unit)
    Integer(I4B), Intent(INOUT) :: unit
    unit = stdout
    Return
  End Subroutine io_getout
  !
  !------------------------------------------------------------------     
  !
  !     Logical unit management
  !
  Subroutine io_assign(lun)
    Integer(I4B), Intent(INOUT) :: lun
    !
    !     Looks for a free unit and assigns it to lun
    !
    Do i = min_lun, max_lun
       If (lun_is_free(i)) Then
          Inquire(unit=i, opened=used, iostat=iostat)
          If (iostat .Ne. 0) used = .True.
          lun_is_free(i) = .False.
          If (.Not. used) Then 
             lun = i
             Return
          Endif
       Endif
    Enddo
    Write(stderr,'(a)') 'No luns available in io_assign'
    Stop 'LUN'
  End Subroutine io_assign
  !
  !===
  !
  Subroutine io_reserve(lun)
    Integer(I4B), Intent(INOUT) :: lun
    !
    !     Useful to specify that one needs to use a particular unit number
    !
    !     For example, assume some legacy code expects to work with unit 15:
    !
    !     call io_reserve(15)   ! this call at the beginning of the program
    !     ...
    !     open(15,....)
    !
    Inquire(unit=lun, opened=used, iostat=iostat)
    If (iostat .Ne. 0) used = .True.
    If (used) Then
       Write(stderr,'(a,i3,a)') &
            'Cannot reserve unit',lun,'. Already connected'
       Stop 'LUN'
    Endif
    If (lun .Ge. min_lun .And. lun .Le. max_lun) lun_is_free(lun) = .False.

    Return
  End Subroutine io_reserve
  !
  !===
  !
  Subroutine io_close(lun)
    Integer(I4B), Intent(INOUT) :: lun
    !
    !     Use this routine instead of a simple close!!
    !
    Close(lun)
    If (lun .Ge. min_lun .And. lun .Le. max_lun) lun_is_free(lun) = .True.
    Return
  End Subroutine io_close
  !
  !===
  !
  Subroutine io_status
    !
    !     Prints a list of the connected logical units and the names of
    !     the associated files
    !

    Write(stdout,'(a)') '******** io_status ********'
    Do i = 0, max_lun
       Inquire(i,opened=opened,named=named,name=filename, &
            form=form,iostat=iostat)
       If (iostat .Eq. 0) Then
          If (opened) Then
             If (named) Then
                Write(stdout,9000) i, form, filename
             Else
                Write(stdout,9000) i, form, 'No name available'
             Endif
          Endif
       Else
          Write(stdout,9000) i, 'Iostat error'
       Endif
    Enddo
    Write(stdout,'(a)') '********           ********'

9000 Format(i4,5x,a,5x,a)
    Return

  End Subroutine io_status
  !
End Module io

!----------------------------------------------------------
!
! FDF (Flexible Data Format) routine package
!
! Copyright Alberto Garcia, Jose Soler, 1996, 1997, 1998
!
!------------------
!
!     Notes: 
!
!     This package implements fdf Version 0.6 (see file fdf.Standard)
!
!     User callable routines:
!
!     fdf_init(filein,fileout)
!     fdf_shutdown
!     fdf_enabled
!     fdf_{integer,single,double,string,boolean} (label,default)
!     fdf_physical(label,default,unit)
!     fdf_block(label,io_unit)
!     fdf_defined(label)
!
!     Implementation notes:
!
!     This version needs the io module for logical unit management.
!
!-------------------------------------------------------------------
!
! Put into F90 module: Chris J. Pickard 1st March 1999
!
!-------------------------------------------------------------------
Module fdf

  Use io

  Implicit None
  !
  Integer, Private, Parameter :: I4B = Selected_int_kind(9)
  Integer, Private, Parameter :: DP  = Kind(1.d0)
  Integer, Private, Parameter :: SP=Kind(1.0)
  !
  Integer(I4B), Private, Parameter :: maxdepth=5
  Integer(I4B), Private :: ndepth, fdf_stack(maxdepth)
  !
  !     Unit numbers for input, output, error notification, and
  !     debugging output (the latter active if fdf_debug is true)
  !
  Integer(I4B), Private :: fdf_in, fdf_out, fdf_err, fdf_log
  Logical, Private      :: fdf_debug, fdf_debug2, fdf_started
  !
  !     Line just read and parsing info
  !
  Character(132), Private :: line
  Integer(I4B), Private, Parameter :: maxntokens = 50
  Integer(I4B), Private :: ntokens
  Integer(I4B), Private :: first(maxntokens), last(maxntokens)
  !
  Data fdf_started, fdf_debug, fdf_debug2 /.False.,.False.,.False./
  !
Contains
  !
  Logical Function fdf_getline()

    Read(fdf_in,End=100,err=100,fmt='(a)') line
    fdf_getline = .True.
    If (fdf_debug2) Write(fdf_log,'(a,a76)') '> ', line
    Call fdf_parse
    Return

100 Continue

    fdf_getline = .False.

    Return
  End Function fdf_getline
  !
  Logical Function leqi(STRNG1,STRNG2)
    !
    !  Case-insensitive lexical equal-to comparison
    !
    Implicit None
    !
    Character*1   :: S1,S2
    Character*(*) :: STRNG1
    Character*(*) :: STRNG2
    !
    Integer(I4B) :: len1, len2, lenc, i
    !
    LEN1=Len(STRNG1)
    LEN2=Len(STRNG2)
    LENC=Min(LEN1,LEN2)
    !
    leqi=.False.
    Do I=1,LENC
       S1=STRNG1(I:I)
       S2=STRNG2(I:I)
       Call CHRCAP(S1,1)
       Call CHRCAP(S2,1)
       If(S1.Ne.S2)Return
    End Do
    ! 
    If(LEN1.Gt.LENC.And.STRNG1(LENC+1:LEN1).Ne.' ')Return
    If(LEN2.Gt.LENC.And.STRNG2(LENC+1:LEN2).Ne.' ')Return
    leqi=.True.
    Return
  End Function leqi
  !  
  Logical Function labeleq(s1,s2)
    !
    !     Compares s1 and s2 without regard for case, or appearance
    !     of '_', '.', '-'.
    !
    Character*(*) :: s1, s2
    Character*80  :: n1, n2

    Call fdf_pack(s1,n1)
    Call fdf_pack(s2,n2)
    labeleq=leqi(n1,n2)
    If (fdf_debug) Then
       If (labeleq .And. .Not. leqi(s1,s2)) Write(fdf_log,'(a,/,a,/,a)') &
            '--------- Considered equivalent:', s1, s2
    Endif
    Return
  End Function labeleq
  !-----------------------------
  !
  Integer(I4B) Function fdf_search(label)
    !
    !     Performs a case-and-punctuation-insensitive search for 'label'
    !     among the tokens in a line.
    !
    !
    Character*(*) :: label
    Integer(I4B) :: i
    !
    fdf_search = 0
    Do i = 1, ntokens
       If (labeleq(label,line(first(i):last(i)))) Then
          fdf_search = i
          Return
       Endif
    Enddo
    !
    Return
  End Function fdf_search
  !-----------------------------------------------------------------
  !
  Logical Function fdf_locate(label)
    !
    !     Searches for label in the fdf hierarchy. If it appears and it
    !     is not part of a comment, the function returns .true. and leaves
    !     the file positioned at the next line. Otherwise, it returns .false.
    !
    !     It supports two kinds of "include" files:
    !
    !     %include filename  
    !     Indicates an unconditional opening of filename for 
    !     further fdf processing.
    !
    !     Label1 Label2 ... < filename  
    !     Indicates that filename should be opened only when 
    !     searching for any of the labels indicated.
    !     'filename' should be an fdf file.
    !

    Character*(*) :: label

    Character*60 :: token1, filename
    Integer(I4B) :: ilabel, iless
    !
    !
    Call fdf_refresh
    If (fdf_debug) Write(fdf_log,'(/,a,1x,a)') 'Looking for ', label

    fdf_locate = .False.

    Rewind(fdf_in)

10  Continue

    If (.Not. fdf_getline()) Then
       If (ndepth .Gt. 1) Then
          Call fdf_close
          Goto 10
       Endif
       If (fdf_debug) Write(fdf_log,'(a,1x,a)') '*Did not find ', label
       Return
    Endif
    !
    If (ntokens .Eq. 0) Goto 10
    !
    token1 = line(first(1):last(1))
    !
    If (leqi(token1,'%include')) Then
       !
       !        Include file
       !
       If (ntokens .Eq. 1) Then
          Write(fdf_err,*) 'FDF: No valid filename after %include'
          Stop
       Endif
       filename = line(first(2):last(2))
       Call fdf_open(filename)
       Goto 10
    Endif

    ilabel = fdf_search(label)
    If (ilabel .Ne. 0) Then
       !
       !        Label found...
       !
       If (leqi(token1,'%block')) Then
          fdf_locate = .True.
          If (fdf_debug) Write(fdf_log,'(a,1x,a)') '*Found ', label
          Return
       Endif

       iless = fdf_search('<')
       If ((iless .Ne. 0) .And. (ntokens .Gt. iless)) Then
          !
          !           Continue search in other file
          !
          filename = line(first(iless+1):last(iless+1))
          Call fdf_open(filename)
          Goto 10
       Endif
       !
       !        If we reach this point we must be dealing with a line
       !        of the form 'Label Value'. But we are not interested if
       !        the string appears in the "Value" section
       !
       If (ilabel .Eq. 1) Then
          fdf_locate = .True.
          If (fdf_debug) Write(fdf_log,'(a,1x,a)') '*Found ', label
          Return
       Else
          Goto 10
       Endif

    Else

       Goto 10

    Endif

  End Function fdf_locate
  !
  Integer(I4B) Function fdf_integer(label,default)
    !
    !     Returns an integer associated with label, or default if label
    !     is not found in the fdf file.
    !

    Character*(*) :: label
    Integer(I4B)  :: default


    !
    Character*10 :: fmtstr
    !
    fdf_integer = default

    If (.Not. fdf_locate(label)) Then
       Write(fdf_out,'(a,5x,i10,5x,a)') label, default, '# Default value'
       Return
    Endif

    If (ntokens.Eq.1) Then
       Write(fdf_err,*) 'FDF_INTEGER: No value for ', label
       Stop
    Endif

    Write(fmtstr,9000) last(2)-first(2)+1
9000 Format('(i',i2.2,')')
    Read(line(first(2):last(2)),fmt=fmtstr) fdf_integer
    Write(fdf_out,'(a,5x,i20)') label, fdf_integer

    Return

  End Function fdf_integer
  !
  Subroutine fdf_init(filein,fileout)
    !
    !     New initialization for fdf. Simplified user interface using
    !     the io package.
    !

    Character*(*) :: filein, fileout

    Integer(I4B) :: debug_level
    !

    !
    !     Prevent the user from opening two head files
    !
    If (fdf_started) Then 
       Write(fdf_err,'(a)') 'FDF: Head file already set...'
       Stop 'HEAD'
    Endif

    Call io_geterr(fdf_err)

    ndepth = 0

    Call fdf_open(filein)
    Call io_assign(fdf_out)
    Open(unit=fdf_out,file=fileout,form='formatted',status='unknown')
    Rewind(fdf_out)

    Write(fdf_out,'(/,a,a,a,i3,/)')&
         '#FDF: Opened ',filein, ' for input. Unit:',fdf_in

    fdf_started = .True.

    debug_level = fdf_integer('fdf-debug',0)
    Call fdf_setdebug(debug_level)

    Return
  End Subroutine fdf_init
  !---------------------------------------------------
  Subroutine fdf_shutdown
    !
    !     Closes the 'head' file
    !


    If (.Not. fdf_started) Return

    Call fdf_refresh
    Call io_close(fdf_in)
    fdf_started = .False.

    Return
  End Subroutine fdf_shutdown
  !--------------------------------------------------
  Logical Function fdf_enabled()

    fdf_enabled = fdf_started
    Return
  End Function fdf_enabled
  !---------------------------------------------------
  Subroutine fdf_setdebug(level)
    !
    !     Debugging levels: 
    !     level <=0: nothing
    !     level  =1: standard
    !     level >=2: exhaustive
    !

    Integer(I4B) :: level
    !
    If (level .Le. 0) Then

       If (fdf_debug) Then
          Call io_close(fdf_log)
          fdf_debug = .False.
       Endif

    Else

       If (.Not. fdf_debug) Then
          Call io_assign(fdf_log)
          Open(fdf_log,file='FDF.debug',form='formatted',status='unknown')
          Rewind(fdf_log)
          fdf_debug = .True.
       Endif
    Endif

    fdf_debug2 = (level .Ge. 2)

    Return
  End Subroutine fdf_setdebug
  !---------------------------------------------------
  !
  Subroutine fdf_open(filename)
    !
    !     Opens a file for fdf processing.
    !
    Character*(*) :: filename


    Integer(I4B) :: lun
    Logical :: file_exists


    ndepth = ndepth + 1
    If (ndepth .Gt. maxdepth) Then
       Write(fdf_err,'(a)') 'FDF: Too many nested fdf files...'
       Stop 'DEPTH'
    Endif

    If (leqi(filename,'stdin')) Then
       lun = 5
       If (fdf_debug) Write(fdf_log,'(a,i1,a)')&
            '--->Reading from Standard Input [depth:', ndepth,'] '

    Else

       Call io_assign(lun)

       Inquire(file=filename,exist=file_exists)
       If (file_exists) Then
          Open(unit=lun,file=filename,status='old',form='formatted')
          Rewind(lun)
          If (fdf_debug) Write(fdf_log,'(a,i1,a,a50)') &
               '--->Opened [depth:', ndepth,'] ', filename
       Else
          Write(fdf_err,'(a,a60)') 'FDF: Cannot open ',filename
       Endif
    Endif

    fdf_stack(ndepth) = lun
    fdf_in = lun

    Return
  End Subroutine fdf_open
  !-----------------------------------
  Subroutine fdf_close
    !
    !     Closes currently opened fdf file, except if it is the original one.
    !
    If (ndepth .Gt. 1) Then
       Call io_close(fdf_in)
       If (fdf_debug) Write(fdf_log,'(a,i1,a)')'--->Closed [depth:', ndepth,']'
       ndepth = ndepth -1
       fdf_in = fdf_stack(ndepth)
    Endif

    Return
  End Subroutine fdf_close
  !-------------------------------------
  Subroutine fdf_refresh
    !
    !     Closes all the open files in the stack (except the first).
    !     Failure to do so would imply that the next Label is searched 
    !     first in the 'deeper' files. fdf_locate calls fdf_refresh 
    !     before doing anything. 
    !

    Integer(I4B) :: i

    Do i = ndepth, 1 , -1
       Call fdf_close
    Enddo

    Return
  End Subroutine fdf_refresh
  !
  Subroutine fdf_parse
    !
    !     Processes the input line looking for meaningful tokens.
    !

    !
    Logical :: intoken, instring

    Integer(I4B) :: c
    Integer(I4B) :: stringdel
    !
    !     Character statement functions
    !
    Integer(I4B) i
    Logical isdigit, isupper, islower, isalpha, isalnum, isextra, istokch
    Logical iscomment, isdelstr, isspecial
    !
    isdigit(i) = (i .Ge. 48) .And. (i .Le. 57)
    isupper(i) = (i .Ge. 65) .And. (i .Le. 90)
    islower(i) = (i .Ge. 97) .And. (i .Le. 122)
    isalpha(i) = isupper(i) .Or. islower(i)
    isalnum(i) = isdigit(i) .Or. isalpha(i)

    !     Extra characters allowed in tokens:  $ % * + & - . / @ ^ _ ~
    isextra(i) = ((i .Ge. 36) .And. (i .Le. 38)) &
         .Or. (i .Eq. 42) .Or. (i .Eq. 43) &
         .Or. (i .Eq. 45) .Or. (i .Eq. 46) &
         .Or. (i .Eq. 47) .Or. (i .Eq. 64) .Or. (i .Eq. 94) &
         .Or. (i .Eq. 95) .Or. (i .Eq. 126)

    istokch(i) = isalnum(i) .Or. isextra(i)
    !
    !     Comments are signaled by:  !  #  ; 
    iscomment(i) = (i.Eq.33) .Or. (i.Eq.35) .Or. (i.Eq.59)
    !
    !     String delimiters: "  '  ` 
    isdelstr(i) = (i.Eq.34) .Or. (i.Eq.39) .Or. (i.Eq.96)
    !
    !     Special characters which are tokens by themselves: <
    isspecial(i) = (i.Eq.60)
    !
    !========================================================
    !
    intoken = .False.
    instring = .False.
    ntokens = 0
    stringdel = 0

    Do i = 1, Len(line)
       c = Ichar(line(i:i))

       If (iscomment(c)) Then
          ! possible comment...
          If (instring) Then
             last(ntokens) = i
          Else
             Goto 1000
          Endif

       Else If (istokch(c)) Then
          ! character allowed in a token...
          If (.Not. intoken) Then
             intoken = .True.
             ntokens = ntokens+1
             first(ntokens) = i
          Endif
          last(ntokens) = i

       Else If (isspecial(c)) Then
          ! character that forms a token by itself...
          If (.Not. instring) Then
             ntokens=ntokens+1
             first(ntokens) = i
             intoken = .False.
          Endif
          last(ntokens) = i

       Else If (isdelstr(c)) Then
          ! string delimiter... make sure it is the right one before
          ! closing the string.
          ! If we are currently in a token, the delimiter is appended to it.

          If (instring) Then
             If (c.Eq.stringdel) Then
                instring = .False.
                intoken = .False.
                stringdel = 0
             Else
                last(ntokens) = i
             Endif
          Else
             If (intoken) Then
                last(ntokens) = i
             Else
                instring = .True.
                stringdel = c
                intoken = .True.
                ntokens = ntokens+1
                first(ntokens) = i+1
                last(ntokens) = i+1
             Endif
          Endif

       Else
          ! token delimiter...

          If (instring) Then
             last(ntokens) = i
          Else
             If (intoken) intoken=.False.
          Endif
       Endif

    Enddo

1000 Continue

    If (fdf_debug2) Then
       Write(fdf_log,*) '            ',  ntokens, ' tokens:'
       Do i=1,ntokens
          Write(fdf_log,*) '                 ','|',line(first(i):last(i)),'|'
       Enddo
    Endif

    Return
  End Subroutine fdf_parse
  !--------------------------------------------------------------------
  !
  !
  !====================================================================
  !
  ! was : character *(*) cjp 1/3/99
  !      character (80) function fdf_string(label,default)
  Character (80) Function fdf_string(label,default)
    !
    !     Returns a string associated with label label, or default if label
    !     is not found in the fdf file.
    !
    Implicit None

    Character*(*) :: label
    Character*(*) :: default

    !
    fdf_string = default

    If (.Not. fdf_locate(label)) Then
       Write(fdf_out,'(a,5x,a,5x,a)') label, default,'# Default value'
       Return
    Endif
    !
    !     From the second token up...
    !
    fdf_string = line(first(2):last(ntokens))
    Write(fdf_out,'(a,5x,a)') label, fdf_string

    Return

  End Function fdf_string
  !
  Logical Function fdf_boolean(label,default)
    !
    !     Returns true if label appears by itself or in the form
    !     label {Yes,true,.true.,T} (case insensitive).
    !
    !     Returns false if label appears in the form
    !     label {No,false,.false.,F} (case insensitive).
    !
    !     If label is not found in the fdf file, fdf_boolean returns the 
    !     logical variable default.
    !

    Character*(*) :: label
    Logical :: default

    Character :: valstr*40

    !
    fdf_boolean = default

    If (.Not. fdf_locate(label)) Then
       Write(fdf_out,'(a,5x,l10,5x,a)') label, default, '# Default value'
       Return
    Endif

    !
    !     If the label appears by itself, we interpret it as .true.
    !
    If (ntokens .Eq. 1) Then
       fdf_boolean = .True.
       Write(fdf_out,'(a,5x,l10,5x,a)') label, fdf_boolean,&
            '# Label by itself'
       Return
    Endif
    !
    !     Look for second word
    !
    valstr=line(first(2):last(2))
    !
    If (leqi(valstr,'yes') .Or. &
         leqi(valstr,'true') .Or. &
         leqi(valstr,'.true.') .Or. &
         leqi(valstr,'t') .Or. &
         leqi(valstr,'y'))       Then

       fdf_boolean = .True.
       Write(fdf_out,'(a,5x,l10)') label, fdf_boolean

    Else If (leqi(valstr,'no') .Or. &
         leqi(valstr,'false') .Or. &
         leqi(valstr,'.false.') .Or. &
         leqi(valstr,'f') .Or. &
         leqi(valstr,'n'))       Then

       fdf_boolean = .False.
       Write(fdf_out,'(a,5x,l10)') label, fdf_boolean

    Else

       Write(fdf_err,*) 'FDF_BOOLEAN: Unexpected fdf logical value ', &
            label, ' = ', valstr
       Stop 

    Endif

    Return
  End Function fdf_boolean
  !
  !----------------------------------------------------------------------
  !

  Real Function fdf_single(label,default)
    !
    !     Returns a single precision value associated with label label, 
    !     or default if label is not found in the fdf file.
    !

    Character*(*) :: label
    Real :: default

    !
    Character*10 :: fmtstr
    !
    fdf_single = default

    If (.Not. fdf_locate(label)) Then
       Write(fdf_out,'(a,5x,g20.10,5x,a)') label, default, '# Default value'
       Return
    Endif

    If (ntokens.Eq.1) Then
       Write(fdf_err,*) 'FDF_SINGLE: No value for ', label
       Stop
    Endif
    Write(fmtstr,9000) last(2)-first(2)+1
9000 Format('(g',i2.2,'.0)')
    Read(line(first(2):last(2)),fmt=fmtstr) fdf_single
    Write(fdf_out,'(a,5x,g20.10)') label, fdf_single

    Return

  End Function fdf_single
  !
  Real(DP) Function fdf_double(label,default)
    !
    !     Returns a double precision value associated with label label, 
    !     or default if label is not found in the fdf file.
    !

    Character*(*) :: label
    Real(DP) :: default



    !
    Character*10 :: fmtstr
    !
    fdf_double = default

    If (.Not. fdf_locate(label)) Then
       Write(fdf_out,'(a,5x,g20.10,5x,a)') label, default, '# Default value'
       Return
    Endif

    If (ntokens.Eq.1) Then
       Write(fdf_err,*) 'FDF_DOUBLE: No value for ', label
       Stop
    Endif
    Write(fmtstr,9000) last(2)-first(2)+1
9000 Format('(g',i2.2,'.0)')
    Read(line(first(2):last(2)),fmt=fmtstr) fdf_double
    Write(fdf_out,'(a,5x,g20.10)') label, fdf_double

    Return

  End Function fdf_double
  !
  !------------------------------------------------------
  Real(DP) Function fdf_physical(label,default,defunit)
    !
    !     Returns a double precision value associated with label label, 
    !     or default if label is not found in the fdf file. Converts
    !     the units to defunit.
    !

    Character*(*) :: label, defunit
    Real(DP) :: default

    Character :: unitstr*10
    Real(DP) :: value

    !
    Character*10 :: fmtstr
    !
    fdf_physical = default

    If (.Not. fdf_locate(label)) Then
       Write(fdf_out,'(a,5x,g20.10,1x,a,5x,a)') &
            label, default, defunit, '# Default value'
       Return
    Endif

    If (ntokens.Eq.1) Then
       Write(fdf_err,*) 'FDF_PHYSICAL: No value for ', label
       Stop
    Endif
    Write(fmtstr,9000) last(2)-first(2)+1
9000 Format('(g',i2.2,'.0)')
    Read(line(first(2):last(2)),fmt=fmtstr) value
    fdf_physical = value
    !
    !     Look for unit
    !
    If (ntokens.Eq.2) Then
       Write(fdf_err,*) 'FDF_PHYSICAL: No unit specified for ', label
       Stop
    Endif

    unitstr=line(first(3):last(3))
    If (.Not. leqi(unitstr,defunit)) &
         fdf_physical = value * fdf_convfac(unitstr,defunit)
    Write(fdf_out,'(a,5x,g20.10,1x,a10)') &
         label, fdf_physical, defunit
    Write(fdf_out,'(a,a,5x,g20.10,1x,a10)') &
         '# Above item originally: ', label, value, unitstr

    Return

  End Function fdf_physical
  !
  !-----------------------------------------------------------------------
  !
  Logical Function fdf_block(label,unit)
    !
    !     Returns "true" and the unit number of the file from which to read
    !     the contents of a block if "label" is associated with a block, and
    !     false if not (unit is set to -1 in this case).
    !
    Character*(*) :: label
    Integer(I4B) :: unit

    Character*50 :: token1, filename
    Integer(I4B) :: iless

    fdf_block = .False.
    unit = -1
    If (.Not. fdf_locate(label)) Return

    token1 = line(first(1):last(1))
    If (.Not. leqi(token1,'%block')) Then
       Write(fdf_err,*) 'FDF_BLOCK: Not a block:',label
       !
       !        Return instead of stopping
       !
       Return
    Endif

    iless = fdf_search('<')
    If ((iless .Ne. 0) .And. (ntokens .Gt. iless)) Then
       !
       !           Read block from file
       !
       filename = line(first(iless+1):last(iless+1))
       If (fdf_debug) Write(fdf_log,'(2a)') '*Reading block from file', &
            filename
       Call fdf_open(filename)
       If (fdf_search('%dump') .Ne. 0) Call fdf_dumpfile(label)
       fdf_block = .True.
       unit = fdf_in
       Return
    Endif
    !
    !     Standard block in fdf file. Dump contents
    !
    Call fdf_dumpblock(label)
    fdf_block = .True.
    unit = fdf_in

    Return
  End Function fdf_block
  !
  !-----------------------------------------------------------------------
  !

  !---------------------------------------
  !     
  Logical Function fdf_defined(label)
    !
    Implicit None

    Character*(*) :: label



    fdf_defined = fdf_locate(label)
    If (fdf_defined) Write(fdf_out,'(a)') label

    Return
  End Function fdf_defined
  !

  !
  !------------------------------------------------------
  !
  Real(DP) Function fdf_convfac( FROM, TO )

    ! Returns conversion factor between a subset of physical units
    ! Written by J.M.Soler. Dec'96.
    ! Modified by Alberto Garcia, Jan'97.

    Character*(*) :: FROM, TO
    Integer(I4B)  :: IU, IFROM, ITO, NU
    Parameter ( NU = 55 )

    Character ::  Dim(NU)*10, NAME(NU)*10
    Real(DP)  :: UNIT(NU)
    !
    !
    !     We allow case variations in the units. This could be dangerous
    !     (meV --> MeV !!) in real life, but not in this restricted 
    !     field.
    !     
    Data (Dim(IU), NAME(IU), UNIT(IU), IU=1,10) / &
         'mass  ', 'Kg      ', 1.D0, &
         'mass  ', 'g       ', 1.D-3, &
         'mass  ', 'amu     ', 1.66054D-27, &
         'length', 'm       ', 1.D0, &
         'length', 'nm      ', 1.D-9, &
         'length', 'Ang     ', 1.D-10, &
         'length', 'Bohr    ', 0.529177D-10, &
         'time  ', 's       ', 1.D0, &
         'time  ', 'fs      ', 1.D-15, &
         'energy', 'J       ', 1.D0/ 
    Data (Dim(IU), NAME(IU), UNIT(IU), IU=11,20) / &
         'energy', 'erg     ', 1.D-7, &
         'energy', 'eV      ', 1.60219D-19, &
         'energy', 'meV     ', 1.60219D-22, &
         'energy', 'Ry      ', 2.17991D-18, &
         'energy', 'mRy     ', 2.17991D-21, &
         'energy', 'Hartree ', 4.35982D-18, &
         'energy', 'K       ', 1.38066D-23, &
         'energy', 'kcal/mol', 6.94780D-21, &
         'force ', 'N       ', 1.D0, &
         'force ', 'eV/Ang  ', 1.60219D-9/
    Data (Dim(IU), NAME(IU), UNIT(IU), IU=21,30) / &
         'force ', 'Ry/Bohr ', 4.11943D-8, &
         'length  ', 'cm      ', 1.d-2, &
         'time    ', 'ps      ', 1.d-12, &
         'time    ', 'ns      ', 1.d-9, &
         'energy  ', 'mHartree', 4.35982D-21, &
         'energy  ', 'kJ/mol  ', 1.6606d-21, &
         'energy  ', 'Hz      ', 6.6262d-34, &
         'energy  ', 'THz     ', 6.6262d-22, &
         'energy  ', 'cm-1    ', 1.986d-23, &
         'energy  ', 'cm^-1   ', 1.986d-23/ 
    Data (Dim(IU), NAME(IU), UNIT(IU), IU=31,40) / &
         'pressure', 'Pa      ', 1.d0, &
         'pressure', 'MPa     ', 1.d6, &
         'pressure', 'GPa     ', 1.d9, &
         'pressure', 'atm     ', 1.01325d5, &
         'pressure', 'bar     ', 1.d5, &
         'pressure', 'Mbar    ', 1.d11, &
         'charge  ', 'C       ', 1.d0, &
         'charge  ', 'e       ', 1.602177d-19, &
         'dipole  ', 'C*m     ', 1.d0, &
         'dipole  ', 'D       ', 3.33564d-30/ 
    Data (Dim(IU), NAME(IU), UNIT(IU), IU=41,50) / &
         'dipole  ', 'debye   ', 3.33564d-30, &
         'dipole  ', 'e*Bohr  ', 8.47835d-30, &
         'dipole  ', 'e*Ang   ', 1.602177d-29,&
         'energy  ', 'cm**-1    ', 1.986d-23, &
         'pressure', 'Ry/Bohr**3', 1.47108d13,&
         'pressure', 'eV/Ang**3 ', 1.60219d11,&
         'MomInert', 'Kg*m**2   ', 1.d0, &
         'MomInert', 'Ry*fs**2  ', 2.17991d-48, &
         'Efield  ', 'V/m       ', 1.d0, &
         'Efield  ', 'V/nm      ', 1.d9 /
    Data (Dim(IU), NAME(IU), UNIT(IU), IU=51,55) / &
         'Efield  ', 'V/Ang     ', 1.d10, &
         'Efield  ', 'V/Bohr    ', 1.8897268d10, &
         'Efield  ', 'Ry/Bohr/e ', 2.5711273d11, &
         'Efield  ', 'Har/Bohr/e', 5.1422546d11, &
         'time    ', 'ps        ', 1.D-12 /
    !
    IFROM = 0
    ITO   = 0
    Do IU = 1,NU
       If (leqi(NAME(IU),FROM)) IFROM = IU
       If (leqi(NAME(IU),TO))   ITO   = IU
    End Do
    If (IFROM .Eq. 0) Then
       Write(6,*) 'FDF_CONVFAC: Unknown unit = ', FROM
       Stop
    Endif
    If (ITO .Eq. 0) Then
       Write(6,*) 'FDF_CONVFAC: Unknown unit = ', TO
       Stop
    Endif

    If (leqi(Dim(IFROM),Dim(ITO))) Then
       FDF_CONVFAC = UNIT(IFROM) / UNIT(ITO)
    Else
       Write(6,*) &
            'FDF_CONVFAC: Unit''s physical dimensions don''t match: ', &
            FROM, ', ', TO
       Stop
    Endif
  End Function fdf_convfac
  !
  !-------------------------------------------------------------------
  !
  Subroutine fdf_dumpblock(label)
    !     
    !     Dumps block contents starting at the current line
    !     

    Character*(*) :: label

    Integer(I4B) :: i, lblock
    Character*60 :: token1

    Write(fdf_out,'(/,a79)') line
    lblock = 0
120 Continue
    If (fdf_getline()) Then
       lblock = lblock + 1
       Write(fdf_out,'(a79)') line
       token1 = line(first(1):last(1))
       If (.Not. leqi(token1,'%endblock')) Goto 120
    Else
       Write(fdf_err,'(a,a,a)') 'FDF_LOCATE: Block ', label, ' does not end!'
       Stop 'FDF'
    Endif
    Write(fdf_out,*)
    !     
    !     Sanity check (optional construct %endblock [ Label [Label] ])
    !     
    If ((ntokens .Gt. 1) .And. &
         (fdf_search(label) .Eq. 0)) Then
       Write(fdf_err,'(a,a,a)') &
            'FDF_LOCATE: Block ', label, ' does not end!'
       Stop 'FDF'
    Endif
    !
    !     Backspace the lines read
    !
    Do i=1,lblock
       Backspace(fdf_in)
    Enddo

    Return
  End Subroutine fdf_dumpblock
  !     
  !--------------------------------------------------------------
  !     
  Subroutine fdf_dumpfile(label)
    !     
    !     Dumps the contents of a file to fdf_out.
    !     The lines are embedded in a %block ... %endblock pair.
    !     


    Character*(*) :: label
    Character :: form*30
    Integer(I4B) :: length


    !     
    !     Build the right format
    !     
    Call chrlen(label,0,length)
    Write(form,'(a,i2.2,a)') '(a,a',length,',10x,a)'

    Write(fdf_out,*)
    Write(fdf_out,form) '%block ', label,'# Originally in include file' 
    !     
    Rewind(fdf_in)
10  Continue
    If (fdf_getline()) Then
       Write(fdf_out,'(a79)') line
       Goto 10
    Endif

    Write(fdf_out,form) '%endblock ', label,'# Originally in include file' 
    Write(fdf_out,*)
    Rewind(fdf_in)
    !     
    Return
  End Subroutine fdf_dumpfile
  !
  !----------------------------------------------------------------
  !
  Subroutine fdf_pack(s,n)
    Implicit None
    Character*(*) :: s, n
    !
    !     Removes occurrences of '_ .-'  from s1
    !
    Character*1 :: c
    Integer(I4B) :: i, j
    Logical :: issep
    issep(i) = (i.Eq.95) .Or. (i.Eq.46) .Or. (i.Eq.45)

    n = ' '
    j = 0
    Do i = 1, Len(s)
       c = s(i:i)
       If (.Not.issep(Ichar(c))) Then
          j = j+1
          n(j:j) = c
       Endif
    Enddo
    Return
  End Subroutine fdf_pack
  !-------------
  Subroutine chrlen(string,nchar,lchar)
    !
    !  CHRLEN accepts a STRING of NCHAR characters and returns LCHAR,
    !  the length of the string up to the last nonblank, nonnull.
    !     
    Character :: char*1
    Character :: string*(*)
    Integer(I4B) :: nchar,lchar
    !
    Integer(I4B) :: ncopy, i
    !
    ncopy=nchar
    If(ncopy.Le.0) ncopy=Len(string)
    !
    Do I=1,ncopy
       lchar=ncopy+1-i
       If(string(lchar:lchar).Ne.' '.And.string(lchar:lchar).Ne.Char(0))Return
    End Do
    lchar=0
    Return
  End Subroutine chrlen
  !
  Subroutine chrcap(string,nchar)
    !
    !  CHRCAP accepts a STRING of NCHAR characters and replaces
    !  any lowercase letters by uppercase ones.
    !
    Integer(I4B) :: nchar, ncopy, i, itemp
    !
    Character :: char*1
    Character :: string*(*)
    !
    ncopy=nchar
    If(ncopy.Le.0) ncopy=Len(string)
    Do I=1,ncopy
       !
       If(Lge(string(i:i),'a').And.Lle(string(i:i),'z'))Then
          itemp=Ichar(string(I:I))+Ichar('A')-Ichar('a')
          string(i:i)=Char(itemp)
       Endif
    End Do
    Return
  End Subroutine chrcap
  !
End Module fdf

