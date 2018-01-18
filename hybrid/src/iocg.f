c $Id: iocg.f,v 1.2 1999/11/26 18:28:14 wdpgaara Exp $

      subroutine iocg( task, naux, cgaux, cgcntr, relaxd, found )

c**************************************************************************
c Reads/writes auxuliary arrays for conj. grag. continuation
c Written by E. Artacho. January 1999.
c**************** INPUT ***************************************************
c character task*3 : 'read' or 'write'
c integer   naux   : dimension of cgaux
c**************** INPUT OR OUTPUT (depending on task) *********************
c real*8    cgaux(naux)   : auxiliary array for CG
c real*8    cgcntr(0:20)  : same
c logical   relaxd        : whether system is relaxed or not.
c***************** OUTPUT *************************************************
c logical found : Has DM been found in disk? (Only when task='read')
c**************************************************************************

      use ionew
      use fdf

      implicit          none
      character         task*(*), paste*33
      logical           found, relaxd
      integer           naux
      double precision  cgaux(naux), cgcntr(0:20)

      external          chkdim, paste


c Internal variables and arrays ------------------------------------------

      character fname*33, sname*30
      logical   exist1, frstme
      integer   nauxr, i, unit1

      save      frstme, fname
      data      frstme /.true./

c ------------------------------------------------------------------------

c find file name ---------------------------------------------------------

      if (frstme) then
        sname = fdf_string('SystemLabel','siesta')
        fname = paste(sname,'.CG')
        frstme = .false.
      endif

c read it if it is there -------------------------------------------------

      if (task.eq.'read' .or. task.eq.'READ') then
        inquire (file=fname, exist=exist1)

        if (exist1) then
          write(6,'(/,a)') 'iocg: Reading CG continuation file'
          call io_assign(unit1)
          open( unit1, file=fname,
     .          form='unformatted', status='unknown' )
          rewind(unit1)
          read(unit1) nauxr, relaxd
          call chkdim( 'iocg', 'cgaux', naux, nauxr, 1 )
          read(unit1) (cgcntr(i), i = 0, 20)
          read(unit1) (cgaux(i), i = 1, naux)
          call io_close(unit1)
          found = .true.
        else
          relaxd = .false.
          cgcntr(0) = 0
          cgcntr(1) = 1
          found = .false.
        endif

c write it ---------------------------------------------------------------

      elseif (task.eq.'write' .or. task.eq.'WRITE') then

        call io_assign(unit1)
        open( unit1, file=fname,
     .        form='unformatted', status='unknown' )
        rewind(unit1)
        write(unit1) naux, relaxd
        write(unit1) (cgcntr(i), i = 0, 20)
        write(unit1) (cgaux(i), i = 1, naux)
        call io_close(unit1)

      else
        stop 'iocg: Incorrect task'
      endif

      return
      end

