c $Id: iofa.f,v 1.2 1999/11/26 18:28:16 wdpgaara Exp $

      subroutine iofa( na, fa )

c *******************************************************************
c Writes forces in eV/Ang
c Emilio Artacho, Feb. 1999
c ********** INPUT **************************************************
c integer na           : Number atoms
c real*8  fa(3,na)     : Forces on the atoms
c *******************************************************************
      use ionew
      use fdf

      implicit          none
      character         paste*33
      integer           na
      double precision  fa(3,*)
      external          paste

c Internal 
      character         sname*30, fname*33
      integer           ia, iu, ix
      logical           frstme
      double precision  Ang, eV
      save              frstme, fname, eV, Ang
      data frstme        /.true./
c -------------------------------------------------------------------

      if (frstme) then
        Ang    = 1.d0 / 0.529177d0
c        eV     = 1.d0 / 13.60580d0
        eV     = 1.d0 / 27.211396132d0
        sname  = fdf_string( 'SystemLabel', 'siesta' )
        fname  = paste( sname, '.FA' )
        frstme = .false.
      endif

c      call io_assign( iu )
c      open( iu, file=fname, form='formatted', status='unknown' )      
c      write(iu,'(i6)') na
c      write(iu,'(i6,3f12.6)') (ia, (fa(ix,ia)*Ang/eV,ix=1,3), ia=1,na)
c      call io_close( iu )
c cambio salida, Nick


      open( 1534, file=fname, form='formatted', status='unknown' )      
      write(1534,'(i6)') na
      write(1534,'(i6,3f18.6)') (ia, (fa(ix,ia)*Ang/eV,ix=1,3), ia=1,na)
      close(1534)

      return
      end
