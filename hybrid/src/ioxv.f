      subroutine ioxv( task, n, u, r, v, fxv ,fv)
c notas, Nick
c task read or write
c n = natot
c u = ucell
c r = rclas
c v = vat
c fxv = foundxv
c fv = foundvat 

c *******************************************************************
c Saves positions and velocities.
c ********** INPUT **************************************************
C  Modules
C
      use precision
      use ionew
      use fdf
      implicit          none

      character         task*(*), paste*33
      logical           fxv,fv
      integer           n
      double precision  r(3,n), v(3,n), u(3,3) 
      external          paste

c Internal variables and arrays
      character  sname*30, fname*33
      integer    ia,iu,iv,ix,na
      logical    frstme
      save       frstme, fname
      data       frstme /.true./

c Find name of file
        if (frstme) then
          sname = fdf_string( 'SystemLabel', 'siesta' )
          fname = paste( sname, '.XV' )
          frstme = .false.
        endif

c Choose between read or write
        if (task.eq.'read' .or. task.eq.'READ') then

c       Check if input file exists
        fxv = .false.

          inquire( file=fname, exist=fxv )
          if (fxv) then

c         Open file
            call io_assign( iu )
            open( iu, file=fname, status='old' )      

c         Read data
             write(6,'(/,a)') 
     .      'ioxv: Reading coordinates and velocities from file'
            do iv = 1,3
              read(iu,*,err=1000,end=1000) (u(ix,iv),ix=1,3)
            enddo
            read(iu,*,err=1000,end=1000) na
            if (na.ne.n) stop 'ioxv: Wrong number of atoms!!'
            do ia = 1,n
              read(iu,*,err=1000,end=1000)
     .          (r(ix,ia),ix=1,3),(v(ix,ia),ix=1,3)
            enddo

c         Close file
            call io_close( iu )

c         Check for velocities
            do ix=1,3
              do ia=1,n
                if(v(ix,ia).ne.0.0) fv=.true.
              enddo
            enddo

          else
c         If input file not found, go to exit point
          write(6,'(/,a)') 'ioxv: WARNING: XV file not found'
            goto 999
          endif

        elseif (task.eq.'write' .or. task.eq.'WRITE') then

c       Open file
          call io_assign( iu )
          open( iu, file=fname, form='formatted', status='unknown' )

c       Write data on file
          write(iu,'(3x,3f18.9)') ((u(ix,iv),ix=1,3),iv=1,3) 
          write(iu,*) n
          do ia = 1,n
            write(iu,'(3f18.9,3x,3f18.9)')
     .        (r(ix,ia),ix=1,3),(v(ix,ia),ix=1,3)
          enddo

c       Close file
          call io_close( iu )

        endif

  999 continue

      return

 1000 stop 'ioxv: problem reading from file'

      end
c====================================================================
      subroutine ioxvconstr( n, u, r, v, s )

c *******************************************************************
c SaveS positions and velocities at each constr step.
c ********** INPUT **************************************************
C  Modules
C
      use precision
      use ionew
      use fdf
      implicit          none

      character         paste*33
      logical           f
      integer           n,s
      double precision  r(3,n), v(3,n), u(3,3) 
      external          paste  

c Internal variables and arrays
      character  sname*30, fname*33,count(1:100)*2
      integer    ia,iu,iv,ix
      data  count/ 
     .     '1','2','3','4','5','6','7','8','9','10',
     .     '11','12','13','14','15','16','17','18','19','20',
     .     '21','22','23','24','25','26','27','28','29','30',
     .     '31','32','33','34','35','36','37','38','39','40',
     .     '41','42','43','44','45','46','47','48','49','50',
     .     '51','52','53','54','55','56','57','58','59','60',
     .     '61','62','63','64','65','66','67','68','69','70',
     .     '71','72','73','74','75','76','77','78','79','80',
     .     '81','82','83','84','85','86','87','88','89','90',
     .     '91','92','93','94','95','96','97','98','99','0'/

      if(s.gt.100) stop 'ioxvcosntr: Increase count variable'
c Find name of file
          sname = fdf_string( 'SystemLabel', 'siesta' )
          fname = paste( sname, '.XV.' )
          fname = paste( fname, count(s) )

c       Open file
          call io_assign( iu )
          open( iu, file=fname, form='formatted', status='unknown' )

c       Write data on file
          write(iu,'(3x,3f18.9)') ((u(ix,iv),ix=1,3),iv=1,3) 
          write(iu,*) n
          do ia = 1,n
            write(iu,'(3f18.9,3x,3f18.9)')
     .        (r(ix,ia),ix=1,3),(v(ix,ia),ix=1,3)
          enddo

c       Close file
          call io_close( iu )

      return
      end

