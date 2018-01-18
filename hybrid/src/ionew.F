      module ionew

c
c Copyright Alberto Garcia, 1996, 1997, 1998, 2000
c
c This module implements an interface to the FORTRAN logical unit
c system. Based on code by Richard Maine.
c
c
c     Logical unit management. Units 0 to min_lun-1 are "reserved",
c     since most of the "typical" files (output, etc) use them.
c
c     Logical units min_lun to min_max are managed by this module.

      use sys
      implicit none
c

      logical, save ::  IOnode

      integer, save ::  stdout, stderr

      integer min_lun, max_lun, nunits
      parameter (min_lun=10, max_lun=99, nunits=max_lun-min_lun+1)
      logical, save ::  lun_is_free(min_lun:max_lun)

      data lun_is_free /nunits*.true./
      data stdout, stderr /6,0/

      CONTAINS

      subroutine io_setup
!
!     Sets IOnode
      ionode = .true.
      end subroutine io_setup

      subroutine io_seterr(unit)
      integer, intent(in) :: unit
      stderr = unit
      end subroutine io_seterr

      subroutine io_setout(unit)
      integer, intent(in) :: unit
      stdout = unit
      end subroutine io_setout


      subroutine io_geterr(unit)
      integer, intent(out) :: unit
      unit = stderr
      end subroutine io_geterr

      subroutine io_getout(unit)
      integer, intent(out) :: unit
      unit = stdout
      end subroutine io_getout

c
c------------------------------------------------------------------     
c
c     Logical unit management
c
      subroutine io_assign(lun)
      integer, intent(out) :: lun

      logical used
      integer iostat

      if (.not. ionode) return
c
c     Looks for a free unit and assigns it to lun
c
      do lun= min_lun, max_lun
         if (lun_is_free(lun)) then
            inquire(unit=lun, opened=used, iostat=iostat)
            if (iostat .ne. 0) used = .true.
            lun_is_free(lun) = .false.
            if (.not. used) return
         endif
      enddo
      write(stderr,'(a)') 'No luns available in io_assign'
      call die('No logical units available')

      end subroutine io_assign
c
      subroutine io_reserve(lun)
      integer, intent(in) :: lun
c
      logical used
      integer iostat

      if (.not. ionode) return

c     Useful to specify that one needs to use a particular unit number
c
c     For example, assume some legacy code expects to work with unit 15:
c
c     call io_reserve(15)   ! this call at the beginning of the program
c     ...
c     open(15,....)
c
      inquire(unit=lun, opened=used, iostat=iostat)
      if (iostat .ne. 0) used = .true.
      if (used) then
             write(stderr,'(a,i3,a)')
     $        'Cannot reserve unit',lun,'. Already connected'
         call die('Cannot reserve unit')
      endif
      if (lun .ge. min_lun .and. lun .le. max_lun)
     $                      lun_is_free(lun) = .false.

      end subroutine io_reserve
c
      subroutine io_close(lun)
      integer, intent(in) :: lun
c
c     Use this routine instead of a simple close!!
c
      close(lun)
      if (lun .ge. min_lun .and. lun .le. max_lun)
     $                     lun_is_free(lun) = .true.

      end subroutine io_close
c
c===
c
      subroutine io_status
c
c     Prints a list of the connected logical units and the names of
c     the associated files
c
      logical named, opened
      integer iostat, i
      character filename*50, form*11

      if (.not. ionode) return

      write(stdout,'(a)') '******** io_status ********'
      do i = 0, max_lun
         inquire(i,opened=opened,named=named,name=filename,
     $           form=form,iostat=iostat)
         if (iostat .eq. 0) then
            if (opened) then
               if (named) then
                  write(stdout,9000) i, form, filename
               else
                  write(stdout,9000) i, form, 'No name available'
               endif
            endif
         else
            write(stdout,9000) i, 'Iostat error'
         endif
      enddo
      write(stdout,'(a)') '********           ********'

 9000 format(i4,5x,a,5x,a)

      end subroutine io_status

      end module ionew






