      module fdf3

      private 
      public fdf_init
      public fdf_get, fdf_block, fdf_bline
      public block, destroy
      public print_block, backspace

      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,100)

c     Declarations for fdf procedures

      interface fdf_get
        module procedure fdf_int, fdf_dp, fdf_bool,
     $                   fdf_sp, fdf_str, fdf_phys
      end interface

      
      interface

         function fdf_defined(label)
         logical fdf_defined
         character(len=*), intent(in) :: label
         end function fdf_defined

         function fdf_enabled()
         logical fdf_enabled
         end function fdf_enabled

         function fdf_block(label,unit)
         logical fdf_block
         character(len=*), intent(in) :: label
         integer, intent(out)  :: unit
         end function fdf_block

         function fdf_convfac(unit1,unit2)
         real*8 fdf_convfac
         character(len=*), intent(in) :: unit1, unit2
         end function fdf_convfac

         subroutine fdf_init(filein,fileout)
         character(len=*), intent(in) :: filein, fileout
         end subroutine fdf_init

         subroutine fdf_inhibit
         end subroutine fdf_inhibit

      end interface

      type line_dlist
         character(len=132)                ::  str
         type(line_dlist), pointer         ::  next
         type(line_dlist), pointer         ::  prev
      end type line_dlist

      type block
      private
         type(line_dlist), pointer         ::  mark
         type(line_dlist), pointer         ::  txt
         type(line_dlist), pointer         ::  last
      end type block

      interface destroy
        module procedure destroy_bp
      end interface

      interface fdf_block
          module procedure fdf_blockf, fdf_blockp
      end interface

      contains

      subroutine destroy_bp(bp)
      type(block), pointer       :: bp
      if (associated(bp%txt)) call destroy_dl(bp%txt)
      deallocate(bp)
      end subroutine destroy_bp

      recursive subroutine destroy_dl(dlp)
      type(line_dlist), pointer       :: dlp
      if (associated(dlp%next)) call destroy_dl(dlp%next)
      deallocate(dlp)
      end subroutine destroy_dl

      subroutine backspace(bp)
      type(block), pointer       :: bp

      if (.not. associated(bp%mark)) then
         bp%mark=>bp%last
      else
         if (.not. associated(bp%mark%prev)) then
            bp%mark => bp%txt
         else
            bp%mark => bp%mark%prev
         endif
      endif
      end subroutine backspace
!
!     Get successive lines from block
!
      function fdf_bline(bp,line) result(res)
      logical res
      type(block), pointer       :: bp
      character(len=*), intent(out)  :: line
      
      res = .false.
 10   continue
         if (.not.associated(bp%mark))      return
         line = bp%mark%str
         bp%mark => bp%mark%next
      if (line(1:1).eq."#") goto 10
      if (line.eq." ") goto 10
      res = .true.
      end function fdf_bline
!
!     Print block
!
      subroutine print_block(bp)
      type(block), pointer       :: bp

      type(line_dlist), pointer       :: p
      character(len=70) :: line

      if (.not. associated(bp)) return
      if (.not. associated(bp%txt)) return
      p=>bp%txt
 5    continue
         if (.not.associated(p)) return
         line = p%str
         write(6,'(a70)') line
         p => p%next
         goto 5
      end subroutine print_block
!
!
!
         function fdf_int(label,default)
         integer fdf_int
         character(len=*), intent(in) :: label
         integer, intent(in) ::  default
         integer fdf_integer
         external fdf_integer
         fdf_int = fdf_integer(label,default)
         end function fdf_int

         function fdf_dp(label,default)
         real(dp) fdf_dp
         character(len=*), intent(in) :: label
         real(dp), intent(in) ::  default
         real(dp) fdf_double
         external fdf_double
         fdf_dp = fdf_double(label,default)
         end function fdf_dp

         function fdf_sp(label,default)
         real(sp) fdf_sp
         character(len=*), intent(in) :: label
         real(sp), intent(in) ::  default
         real(sp) fdf_single
         external fdf_single
         fdf_sp = fdf_single(label,default)
         end function fdf_sp

         function fdf_phys(label,default,unit)
         real(dp) fdf_phys
         character(len=*), intent(in) :: label, unit
         real(dp), intent(in) ::  default
         real(dp) fdf_physical
         external fdf_physical
         fdf_phys = fdf_physical(label,default,unit)
         end function fdf_phys

         function fdf_bool(label,default)
         logical fdf_bool
         character(len=*), intent(in) :: label
         logical, intent(in) ::  default
         logical fdf_boolean
         external fdf_boolean
         fdf_bool = fdf_boolean(label,default)
         end function fdf_bool

         function fdf_str(label,default)
         character(len=80) fdf_str
         character(len=*), intent(in) :: label
         character(len=*), intent(in) ::  default
         character(len=80) fdf_string
         external fdf_string
         fdf_str =  fdf_string(label,default)
         end function fdf_str

         function fdf_blockf(label,unit)
         logical fdf_blockf
         character(len=*), intent(in) :: label
         integer, intent(out) ::  unit
         logical fdf_block
         external fdf_block
         fdf_blockf = fdf_block(label,unit)
         end function fdf_blockf
!
!        Fill in block structure
!
         function fdf_blockp(label,bp) result(res)
         logical res
         character(len=*), intent(in) :: label
         type(block), pointer         :: bp

         integer unit
         character(len=132) line
         logical head
         type(line_dlist), pointer         :: p

         if (associated(bp)) call destroy(bp)
         head = .true.

         res = fdf_blockf(label,unit)
         if (res) then
            allocate(bp)
            nullify(bp%mark)
            nullify(bp%txt)
 1          continue
            read(unit,fmt='(a)',iostat=ierr) line
            if (ierr .eq. 0) then
               if (line(1:9).ne."%endblock") then
                  if (head) then
                     allocate(bp%txt)
                     nullify(bp%txt%prev)
                     p => bp%txt
                     bp%mark=>bp%txt
                     bp%last=>bp%txt
                     head = .false.
                  else
                     allocate(p%next)
                     p%next%prev => p
                     p=>p%next
                  endif
                  p%str = line
                  nullify(p%next)
                  bp%last => p
                  goto 1
               endif
            endif
            if (.not. associated(bp%txt)) then
              !!! Empty block!!!
              res = .false.
            endif
         endif
         end function fdf_blockp

      end module fdf3





