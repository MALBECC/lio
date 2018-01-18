      SUBROUTINE timestamp (str)

      implicit none

      CHARACTER(len=*) str

      character(len=1), parameter :: dash = "-"
      character(len=1), parameter :: colon = ":"
      character(len=3), parameter :: prefix = ">> "
      character(len=3)  :: month_str(12)
      data month_str 
     $    /'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP',
     &     'OCT','NOV','DEC'/

      integer sec, min, hour, day, month, year
      integer values(8)

      call date_and_time(values=values)
      year = values(1)
      month = values(2)
      day = values(3)
      hour = values(5)
      min = values(6)
      sec = values(7)
      
      write(6,1000) prefix, trim(str), colon,
     $              day, dash, month_str(month), dash, year,
     $              hour, colon, min, colon, sec

 1000 format(2a,a1,2x,i2,a1,a3,a1,i4,2x,i2,a1,i2.2,a1,i2.2)

      END

