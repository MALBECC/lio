subroutine vec_init(Ene,Vec,N,vecnum,Mlr,NCO,Nvirt,Ndim,Sdim,maxIt)
use garcha_mod  , only: npas
use excited_data, only: use_last, guessLR, estda
   implicit none

   integer, intent(in)   :: N, Mlr, NCO, Nvirt, Ndim
   integer, intent(inout):: vecnum, Sdim, maxIt
   LIODBLE, intent(in)   :: Ene(Mlr)
   LIODBLE, intent(out)  :: Vec(N,vecnum)

   logical :: find
   integer :: ii, occ, virt, cont, start
   integer, dimension(:), allocatable :: ind
   LIODBLE, dimension(:), allocatable :: deltaE

   if ( npas > 1 .and. use_last ) then
      print*, "Using last step as initial guess"
      if (.not. allocated(guessLR)) then
         print*, "The guessLR array has been not allocated"
         stop
      endif
      Vec = guessLR
   else
      print*, "Generating Trial Vectors as initial guess"
      allocate(deltaE(Ndim),ind(Ndim))
      ! Calculate delta molecular orbital energies
      cont = 1
      do ii=1,Ndim
         cont = ii - 1
         occ = NCO - (cont/Nvirt)
         virt = mod(cont,Nvirt) + NCO + 1
         deltaE(ii) = Ene(virt) - Ene(occ)
      enddo

      ! Sorting Energies
      ind = 0
      call eigsort(deltaE,ind,Ndim)

      ! Energy-Specific TDA. DOI:10.1021/ct200485x
      start = 1
      if (estda) call estda_window(deltaE,ind,Ndim,start,vecnum,Sdim,maxIt)

      Vec = 0.0D0
      do ii=1,vecnum
         Vec(ind(start),ii) = 1.0D0
         start = start + 1
      enddo
      deallocate(deltaE,ind)
   endif
end subroutine vec_init

subroutine eigsort(a,b,N)

   implicit none

   integer, intent(in) :: N
   LIODBLE, intent(in):: a(N)
   integer, intent(out) :: b(N)
   integer :: i,imin, temp_i
   LIODBLE :: temp_r
   LIODBLE :: a2(N)

   a2 = a
   do i = 1, N
      b(i) = i
   end do
   do i = 1, N-1
      imin = minloc(a2(i:),1) + i - 1
      if (imin /= i) then
         temp_r = a2(i); a2(i) = a2(imin); a2(imin) = temp_r
         temp_i = b(i) ; b(i)  = b(imin) ; b(imin)  = temp_i
      end if
   end do
end subroutine

