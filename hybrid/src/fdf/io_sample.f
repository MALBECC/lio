c
      program io_sample

      integer stderr
      
      call io_geterr(stderr)

      call io_assign(lun)
      open(lun,file='OUTPUT',status='unknown')
      call io_setout(lun)

      call io_status

      call io_assign(lun)
      open(lun,file='file2',status='unknown')

      call io_status

      call io_close(lun)

      call io_assign(lun)
      open(lun,file='file3',form='unformatted',status='unknown')

      call io_status

      write(stderr,*) ' Error messages appear here'

      end


