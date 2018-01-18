      module sys

      CONTAINS

      subroutine die(str)
      character(len=*), intent(in), optional   :: str
      integer :: Node
      Node = 0
      if (Node.eq.0) then
         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
         write(6,'(a)') 'Stopping Program'
      endif
      stop
      end subroutine die

      end module sys

