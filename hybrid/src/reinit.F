      subroutine reinit(slabel, sname) 

C Subroutine to initialise the reading of the data for SIESTA 
C
C     It uses the FDF (Flexible Data Fromat) package 
C     of J.M.Soler and A.Garcia
C
C Taken from redata. Writen by P.Ordejon, December'96
C **************************** OUTPUT *********************************
C character*20 slabel      : System Label (to name output files)
C character*59 sname       : System Name
C **********************************************************************

C
C  Modules
C
      use precision
      use ionew
      use fdf

      implicit none

      character
     .  slabel*20, sname*59, string*20

C  Internal variables .................................................
      character
     .  filein*20, fileout*20, line*150 

      character
     .  slabel_default*59, sname_default*20

      integer 
     .  count, Node, length, lun, lun_tmp

      logical debug_input, file_exists


C Get Node number
      Node = 0

C Print Welcome and Presentation .......................................
      if (Node.eq.0) then
        write(6,'(/a)') 
     .   '                           ***********************       '
        write(6,'(a)') 
     .   '                           *  WELCOME TO HYBRID  *       '
        write(6,'(a)')
     .   '                           ***********************       '
C ..................

C Dump data file to output file 
C and generate scratch file for FDF to read from
C (except if INPUT_DEBUG exists)
C

        inquire(file='INPUT_DEBUG',exist=debug_input)
        if (debug_input) then
           write(6,'(a)') 'WARNING: ' //
     $          'Siesta is reading its input from file INPUT_DEBUG'
           
           call io_assign(lun)
           filein = 'INPUT_DEBUG'
           open(lun,file='INPUT_DEBUG',form='formatted',status='old')
           rewind(lun)
        else

!          Read from standard input

           write(6,'(/a)') 'reinit: Reading from standard input'
           lun = 5
           call io_assign(lun_tmp)
           do
             call system_clock( count )
             write(string,*) count
             filein = 'INPUT_TMP.'//adjustl(string)
             inquire( file=filein, exist=file_exists )
             if (.not.file_exists) exit
           end do
           open(lun_tmp,file=filein,
     $                  form='formatted',status='replace')
           rewind(lun_tmp)
        endif

        write(6,'(a,23(1h*),a,28(1h*))')
     .    '***', ' Dump of input data file '
!
 10     continue
        read(lun,err=20,end=20,fmt='(a)') line
        call chrlen(line,0,length)
        if (length .ne. 0) then
           write(6,'(a)') line(1:length)
           if (.not. debug_input) write(lun_tmp,'(a)') line(1:length)
        endif
        goto 10
 20     continue
!
        write(6,'(a,23(1h*),a,29(1h*))')
     .    '***', ' End of input data file '
!
!       Choose proper file for fdf processing
!
        if (debug_input) then
           call io_close(lun)
        else
           call io_close(lun_tmp)
        endif

! Set up fdf ...
        fileout = 'out.fdf'
        call fdf_init(filein,fileout)

! Defile Name of the system ...
        sname_default = ' '
        sname = fdf_string('SystemName',sname_default)
        write(6,'(/a,71(1h-))') 'reinit: '
        write(6,'(a,a)') 
     .   'reinit: System Name: ',sname
        write(6,'(a,71(1h-))') 'reinit: '
! ...

! Defile System Label (short name to label files) ...
        slabel_default  = 'Scarlett'
        slabel = fdf_string('SystemLabel',slabel_default)
        write(6,'(a,a)') 
     .   'reinit: System Label: ',slabel
        write(6,'(a,71(1h-))') 'reinit: '
! ...
      endif

      return
      end

