!##########################################################
!##########################################################
module trans_Data
   implicit none
   real*8, dimension(:,:), allocatable :: rho_exc
   logical :: gaussian_convert=.false.
   contains

   subroutine translation(M,Mat_lio)

   implicit none
  
   integer, intent(in) :: M
   real*8, dimension(M,M), intent(out) :: Mat_lio

   integer :: MM, i, j, idx, l, h, stot, ptot, dtot
   real*8, dimension(:,:), allocatable :: Matrix, Mat_par
   real*8, dimension(:), allocatable :: Vector
   integer, dimension(:), allocatable :: gaussian_ind, indices
   character(len=4) :: mom

   print*, "" 
   print*, "#####################################" 
   print*, "        SUBROUTINE TRADUCTION        " 
   print*, "#####################################" 
   print*, "" 

! Reading the Gaussian input
   open(unit=457, file="gau_inp")
     allocate (gaussian_ind(M))
      do i=1,M
        read(457,*) mom
        if ((INDEX(mom,"S"))>0) gaussian_ind(i) = 0
        if ((INDEX(mom,"P"))>0) gaussian_ind(i) = 1
        if ((INDEX(mom,"D"))>0) gaussian_ind(i) = 2
      enddo
   close(457)

   MM=M*(M+1)/2

   allocate(Matrix(M,M),Mat_par(M,M),Vector(MM),indices(M))

   open(unit=456,file="gau_dens")
     do i=1, MM
        read(456,*) Vector(i)
     enddo
   close(456)

   idx=0
   do j = 1, M
     do i = 1, j - 1
       idx = i + (j*(j-1)/2)
       Matrix(i,j) = Vector(idx)
       Matrix(j,i) = Vector(idx)
     enddo
       idx = j + (j*(j-1)/2)
       Matrix(j,j) = Vector(idx)
   enddo

   stot=0
   ptot=0
   dtot=0

   do i=1,M
      if (gaussian_ind(i) == 0) then
         stot = stot + 1
      else if (gaussian_ind(i) == 1) then
         ptot = ptot + 1
      else
         dtot = dtot + 1
      endif
   enddo

   j=1
   l=1
   h=1

   do i=1,M
      if(gaussian_ind(i) == 0) then
        indices(j) = i
        j = j + 1
      else if (gaussian_ind(i) == 1) then
        indices(stot + l) = i
        l = l + 1
      else
        indices(ptot + h) = i
        h = h + 1
      endif
   enddo
 
! Obtain the Density matrix in LIO format
   do i=1,M
     Mat_par(i,:) = Matrix(indices(i),:)
   enddo

   Mat_lio = Mat_par

   do j=1,M
      Mat_lio(:,j) = Mat_par(:,indices(j))
   enddo

   end subroutine translation
end module trans_Data
!##########################################################
!##########################################################
