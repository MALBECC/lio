program average_spectra
   implicit none

   character(len=50) :: input1, input2, input3
   integer           :: iStep, nSteps
   real(kind=8)      :: crd
   real(kind=8), allocatable :: mu1(:), mu2(:), mu3(:)


   call get_command_argument(1, input1)
   call get_command_argument(2, input2)
   call get_command_argument(3, input3)

   open(101, file=input1)
   open(102, file=input2)
   open(103, file=input3)
   open(104, file='esp_sum')
          
   call obtain_lines(nSteps, input1)

   allocate(mu1(1:nSteps))
   allocate(mu2(1:nSteps))
   allocate(mu3(1:nSteps))

   do iStep = 1, nSteps
      read(101,*)  crd, mu1(iStep)
      read(102,*)  crd, mu2(iStep)
      read(103,*)  crd, mu3(iStep)
      write(104,*) crd, (mu1(iStep) + mu2(iStep) + mu3(iStep)) / 3
   enddo

   close(101)
   close(102)
   close(103)
   close(104)
   deallocate(mu1, mu2, mu3)


contains

   ! Gets the number of lines in a file.
   subroutine obtain_lines(n_lines, input_file)
      implicit none
      character(len=50), intent(in)    :: input_file
      integer          , intent(inout) :: n_lines
   
      integer            :: ios
   
      n_lines = 0
      open(100, file = input_file)
      do
         read(100,*,iostat = ios)
         if (ios /= 0) exit
         n_lines = n_lines + 1
      enddo
      close(100)
   endsubroutine obtain_lines

end program average_spectra
