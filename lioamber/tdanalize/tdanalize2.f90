
          implicit none

          double precision, allocatable :: mu1(:), mu2(:), mu3(:)
          double precision a
          integer i
          integer steps

          character(30) input1, input2, input3


          call get_command_argument(1,input1)
          call get_command_argument(2,input2)
          call get_command_argument(3,input3)

          open (unit=1,file=input1)
          open (unit=2,file=input2)
          open (unit=3,file=input3)
          open (unit=4, file='esp_suma')
          
          write(*,*) 'steps'
          read(*,*) steps


          allocate (mu1(1:steps))
          allocate (mu2(1:steps))
          allocate (mu3(1:steps))

          do i = 1, steps
           read(1,*) a,mu1(i)
           read(2,*) a,mu2(i)
           read(3,*) a,mu3(i)
           write (4,*) a,(mu1(i)+mu2(i)+mu3(i))/3
          enddo
          deallocate(mu1); deallocate(mu2); deallocate(mu3)
          end
