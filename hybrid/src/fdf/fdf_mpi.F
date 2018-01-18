      module fdf
#ifdef MPI
      use mpi
#endif
!
!  Copyright Alberto Garcia, Jose M. Soler (1996-)
!
!=====================================================================
!
!     This module implements an extended Fortran 90/95 interface
!     to the Flexible Data Format library of A. Garcia and J.M. Soler,
!     originally written in Fortran 77.
!
!     It provides interface blocks for the routines in the f77 library,
!     as well as new routines which take advantage of the new facilities
!     in Fortran 90/95.
!
!
!     NEW FEATURES:
!
!     a) Block pointers. 
!
!     Block content can now be flexibly handled by means of a pointer
!     to a derived type "block". Typical usage:
!
!     use fdf
!     type(block), pointer :: bp
!
!     if (fdf_block('SomeBlock',bp)) then
!         loop: do
!               if (.not. fdf_bline(bp,line)) exit loop
!               (process line, possibly with 'digest')
!         enddo loop
!     endif
!     call destroy(bp)
!
!     The generic name 'fdf_block' stands for both the old and the new
!     syntax (In the old syntax: fdf_block('SomeBlock',unit), the function
!     returns in 'unit' a unit number from which to read the contents of
!     the block.)
!
!     Routine fdf_bline returns in 'line' the next non-blank, 
!     non-comment line from the block, unless there are no more
!     lines, in which case it returns .false. and 'line' is undefined.
!     Optionally, by specifying 'any=.true.', fdf_bline will return
!     the next line, even if blank.
!
!     Routine 'backspace' moves an internal pointer to the previous line
!     returned (unless the optional 'last_physical=.true.' is specified, 
!     in which case the internal pointer will move to the previous physical 
!     line in the block (even if a blank or comment line); this can lead 
!     to unpredictable behavior) 
!
!     Routine 'rewind' moves the internal pointer to the beginning of 
!     the block.
!
!     The call to 'destroy' frees the storage associated to the pointer bp.
!
!     Among the advantages of the new interface to block handling, are:
!
!     * Automatic detection of the end of the block.
!     * Block pointers can be kept around for as long as needed.
!     * Will work even if the underlying mechanism for block management
!       is not based on files.
!
!     b) Generic interface to scalar routines.
!
!     The generic function 'fdf_get' can be used instead of any of the
!     old scalar routines. The old names are also accepted, for backwards
!     compatibility. 
!
!     c) New routine returning a "digested" string.
!
!     Function fdf_parsed_string(label,default) returns a pointer to
!     a 'parsed_line' derived type (see module parse).
!
!------------------------------------------------------------------------
!     Future enhancements:
!
!     * Stand-alone f90/95 version.
!     * Better exception handling.
!     * MPI-awareness.
!     * Array primitives.
!     * Nested FDF namespaces.
!     * Non-recursive implementation of 'destroy'
!========================================================================

      implicit none

      private 

#ifdef MPI
      integer MPIerror
#endif
      integer Node, Nodes
!
!     Routines already in the f77 library
!
      public fdf_init
      public fdf_convfac
      public fdf_integer, fdf_single, fdf_double, fdf_physical
      public fdf_string, fdf_boolean
      public fdf_defined, fdf_enabled
!
      public fdf_block
!
!     New routines
!
      public block, destroy
      public print_block, backspace, rewind
      public fdf_bline
      public fdf_get
      public fdf_parsed_string
!
!     Private kind declarations
!
      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,100)

c     Declarations for MPI-wrapped fdf procedures

      interface fdf_init
        module procedure fdf_init_wrap
      end interface
      interface fdf_integer
        module procedure fdf_integer_wrap
      end interface
      interface fdf_single
        module procedure fdf_single_wrap
      end interface
      interface fdf_double
        module procedure fdf_double_wrap
      end interface
      interface fdf_physical
        module procedure fdf_physical_wrap
      end interface
      interface fdf_boolean
        module procedure fdf_boolean_wrap
      end interface
      interface fdf_string
        module procedure fdf_string_wrap
      end interface
      interface fdf_defined
        module procedure fdf_defined_wrap
      end interface
      interface fdf_enabled
        module procedure fdf_enabled_wrap
      end interface

      interface fdf_get
        module procedure fdf_integer_wrap, fdf_single_wrap,
     $                   fdf_double_wrap, fdf_boolean_wrap,
     $                   fdf_string_wrap, fdf_physical_wrap
      end interface

      interface fdf_block
          module procedure fdf_blockp, fdf_blockf_wrap
      end interface

      interface destroy
        module procedure destroy_bp
      end interface

      interface backspace
        module procedure backspace_fdf_block
      end interface

      interface rewind
        module procedure rewind_fdf_block
      end interface

!     New derived types to support blocks
!
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
         type(line_dlist), pointer         ::  last_line_returned
      end type block

      interface

         function leqi(s1,s2)
         logical leqi
         character(len=*), intent(in)   :: s1, s2
         end function leqi

         function fdf_convfac(unit1,unit2)
         real(selected_real_kind(14,100)) fdf_convfac
         character(len=*), intent(in) :: unit1, unit2
         end function fdf_convfac

      end interface

      CONTAINS

!-------------------------------------------------------------------
      subroutine fdf_init_wrap(filein,fileout)
      character(len=*), intent(in) :: filein, fileout
      external fdf_init
#ifdef MPI
      call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
      call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
#else
      Node = 0
      Nodes = 1
#endif
      if (Node.eq.0) call fdf_init(filein,fileout)
      end subroutine fdf_init_wrap

      function fdf_defined_wrap(label) result(defined)
      logical defined
      character(len=*), intent(in) :: label
      logical fdf_defined
      external fdf_defined
      if (Node.eq.0) defined = fdf_defined(label)
#ifdef MPI
      call MPI_Bcast(defined,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
      end function fdf_defined_wrap

      function fdf_enabled_wrap() result(enabled)
      logical enabled
      logical fdf_enabled
      external fdf_enabled
      if (Node.eq.0) enabled = fdf_enabled()
#ifdef MPI
      call MPI_Bcast(enabled,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
      end function fdf_enabled_wrap

      function fdf_integer_wrap(label,default) result(int)
      integer int
      character(len=*), intent(in) :: label
      integer, intent(in) ::  default
      integer fdf_integer
      external fdf_integer
      if (Node.eq.0) int = fdf_integer(label,default)
#ifdef MPI
      call MPI_Bcast(int,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
      end function fdf_integer_wrap

      function fdf_single_wrap(label,default) result(single)
      real single
      character(len=*), intent(in) :: label
      real, intent(in) ::  default
      real fdf_single
      external fdf_single
      if (Node.eq.0) single = fdf_single(label,default)
#ifdef MPI
      call MPI_Bcast(single,1,DAT_single,0,MPI_Comm_World,MPIerror)
#endif
      end function fdf_single_wrap

      function fdf_double_wrap(label,default) result(double)
      real*8 double
      character(len=*), intent(in) :: label
      real*8, intent(in) ::  default
      real*8 fdf_double
      external fdf_double
      if (Node.eq.0) double = fdf_double(label,default)
#ifdef MPI
      call MPI_Bcast(double,1,DAT_double,0,MPI_Comm_World,MPIerror)
#endif
      end function fdf_double_wrap

      function fdf_physical_wrap(label,default,unit) result(physical)
      real*8 physical
      character(len=*), intent(in) :: label, unit
      real*8, intent(in) ::  default
      real*8 fdf_physical
      external fdf_physical
      if (Node.eq.0) physical = fdf_physical(label,default,unit)
#ifdef MPI
      call MPI_Bcast(physical,1,DAT_double,0,MPI_Comm_World,MPIerror)
#endif
      end function fdf_physical_wrap

      function fdf_boolean_wrap(label,default) result(boolean)
      logical boolean
      character(len=*), intent(in) :: label
      logical, intent(in) ::  default
      logical fdf_boolean
      external fdf_boolean
      if (Node.eq.0) boolean = fdf_boolean(label,default)
#ifdef MPI
      call MPI_Bcast(boolean,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
      end function fdf_boolean_wrap

      function fdf_string_wrap(label,default) result(string)
      character(len=132) string
      character(len=*), intent(in) :: label
      character(len=*), intent(in) ::  default
      character(len=132) fdf_string
      external fdf_string
      if (Node.eq.0) string = fdf_string(label,default)
#ifdef MPI
      call MPI_Bcast(string,132,MPI_character,0,MPI_Comm_World,MPIerror)
#endif
      end function fdf_string_wrap

!!!----------------------------------------------------------------
      subroutine destroy_bp(bp)
      type(block), pointer       :: bp
      type(line_dlist), pointer       :: dlp, dlp2
      if (associated(bp%txt)) call destroy_dl(bp%txt)
      deallocate(bp)
      end subroutine destroy_bp

!-------------------------------------------------------------------
      recursive subroutine destroy_dl(dlp)
      type(line_dlist), pointer       :: dlp
      if (associated(dlp%next)) call destroy_dl(dlp%next)
      deallocate(dlp)
      end subroutine destroy_dl

!-------------------------------------------------------------------
!
      function fdf_parsed_string(label,default) result(pline)
      use parse
      type(parsed_line), pointer ::  pline
      character(len=*), intent(in) :: label
      character(len=*), intent(in) ::  default

      character(len=132) temp
      temp = fdf_string(label,default)
      pline=>digest(temp)
      end function fdf_parsed_string
!
!===================================================================

         function fdf_blockf_wrap(label,unit)
         logical fdf_blockf_wrap
         character(len=*), intent(in) :: label
         integer, intent(out) ::  unit

         interface         ! for the external fdf_block...
            function fdf_block(label,unit)
            logical fdf_block
            character(len=*), intent(in) :: label
            integer, intent(out)  :: unit
            end function fdf_block
         end interface

         if (Node.ne.0)
     $        call die("Cannot use old block interface in MPI")
         fdf_blockf_wrap = fdf_block(label,unit)
         end function fdf_blockf_wrap

         function fdf_blockf(label,unit)
         logical fdf_blockf
         character(len=*), intent(in) :: label
         integer, intent(out) ::  unit

         interface         ! for the external fdf_block...
            function fdf_block(label,unit)
            logical fdf_block
            character(len=*), intent(in) :: label
            integer, intent(out)  :: unit
            end function fdf_block
         end interface

         fdf_blockf = fdf_block(label,unit)
         end function fdf_blockf
!
!-------------------------------------------------------------------
!        Fill in block structure
!
         function fdf_blockp(label,bp) result(res)
         use parse
         logical res
         character(len=*), intent(in) :: label
         type(block), pointer         :: bp

         integer unit, ierr
         character(len=132) line
         logical head
         type(line_dlist), pointer         :: dlp
         type(parsed_line), pointer        :: p

         if (associated(bp)) call destroy(bp)
         head = .true.

         if (Node.eq.0) res = fdf_blockf(label,unit)
#ifdef MPI
         call MPI_Bcast(res,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
         if (.not.res) return

         allocate(bp)
         nullify(bp%mark)
         nullify(bp%txt)

         loop: DO
            if (Node.eq.0)
     $            read(unit,fmt='(a)',iostat=ierr) line
#ifdef MPI
            call MPI_Bcast(ierr,1,MPI_integer,0,
     $                            MPI_Comm_World,MPIerror)
#endif
            if (ierr .ne. 0) exit loop
#ifdef MPI
            call MPI_Bcast(line,132,MPI_character,0,
     $                            MPI_Comm_World,MPIerror)
#endif

           p=>digest(line)
           if (ntokens(p) .ge. 1) then
              if (leqi(tokens(p,1),"%endblock")) exit loop
           endif

           if (head) then
              allocate(bp%txt)
              nullify(bp%txt%prev)
              dlp => bp%txt
              bp%mark=>bp%txt
              bp%last=>bp%txt
              head = .false.
           else
              allocate(dlp%next)
              dlp%next%prev => dlp
              dlp=>dlp%next
           endif

           dlp%str = line
           nullify(dlp%next)
           bp%last => dlp
           call destroy(p)

        enddo loop

        if (.not. associated(bp%txt)) then
           !!! Empty block!!!
           call warn("fdf_blockp: Block is empty...")
           res = .false.
        else
           bp%last_line_returned => bp%txt   
        endif

         end function fdf_blockp
!
!-------------------------------------------------------------------
      subroutine backspace_fdf_block(bp,physical_line)
      type(block), pointer       :: bp
      logical, intent(in), optional :: physical_line

      logical last_physical_line

      last_physical_line = .false.
      if (present(physical_line)) last_physical_line = physical_line

      if (.not. last_physical_line) then
!        Put the mark at the point of the last returned line,
!        as determined in fdf_bline
         bp%mark => bp%last_line_returned
         return
      endif
!
!     Backspace to the previous physical line in the block
!     (i.e., it might be a blank or a comment line)
!
      if (.not. associated(bp%mark)) then  ! We are at the end of block
         bp%mark=>bp%last
      else
         if (.not. associated(bp%mark%prev)) then ! at the beginning
            bp%mark => bp%txt
         else
            bp%mark => bp%mark%prev
         endif
      endif
      end subroutine backspace_fdf_block

!-------------------------------------------------------------------
      subroutine rewind_fdf_block(bp)
      type(block), pointer       :: bp

      if (.not. associated(bp))
     $     call die("rewind: Block not associated")
      if (.not. associated(bp%txt))
     $     call die("rewind: Block text not associated")
      bp%mark=>bp%txt
      bp%last_line_returned=>bp%txt
      end subroutine rewind_fdf_block
!
!
!-------------------------------------------------------------------
!     Get successive non-blank, non-comment lines from block
!     Optionally, if 'any=.true.' is specified, return any line,
!     even if blank or a comment line.
!
      function fdf_bline(bp,line,any) result(res)
      use parse
      logical res
      type(block), pointer       :: bp
      character(len=*), intent(out)  :: line
      logical, intent(in), optional  :: any
      
      type(parsed_line), pointer   :: p
      logical any_line

      any_line = .false.
      if (present(any)) any_line = any

      res = .false.
      loop: do
         if (.not.associated(bp%mark))      return
         line = bp%mark%str
         bp%last_line_returned=>bp%mark
         bp%mark => bp%mark%next
         p => digest(line)
         if ((ntokens(p) .ne. 0)  .or. any_line) exit loop
      enddo loop
      call destroy(p)
      res = .true.
      end function fdf_bline

!-------------------------------------------------------------------
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
!-------------------------------------------------------------------
!
!     Private copy of warn and die
!
      subroutine warn(str)
      character(len=*) :: str
      if (Node.eq.0) write(6,'(2a)') '*WARNING: ', str
      end subroutine warn

      subroutine die(str)
      character(len=*), intent(in), optional   :: str

      if (Node.eq.0) then
         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
         write(6,'(a)') 'Stopping Program'
         write(6,*) 1.0/real(Node)
      endif
#ifdef MPI
      call MPI_Abort(MPI_Comm_World,1,MPIerror)
!!      call MPI_Finalize(rc)
      stop
#else
      stop
#endif
      end subroutine die

      end module fdf





