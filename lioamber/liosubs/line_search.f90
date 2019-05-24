subroutine line_search(n_points, Energy, step_size, lambda)
   ! Performs a linear search for a given Energy array.
   !
   !  If minimum value of E is E(n_points)
   !     return lambda = step_size * n_points
   !  If minimum value of E is E(1)
   !     return lambda = 0.d0
   !  else
   !     returns an expected lambda that minimice E(lambda)
   !     using a parabolic interpolation
   implicit none
   integer     , intent(in)  :: n_points
   real(kind=8), intent(in)  :: Energy(n_points), step_size
   real(kind=8), intent(out) :: lambda

   integer      :: min_Energy_position
   real(kind=8) :: dE1, dE2, modif_fac

   if (n_points <= 2) then
      write(*,*) "Wrong n_points in lineal search, n_points need to be > 2"
      stop
   end if

   ! Finds minimum value in Energy elements
   min_energy_position = minloc(Energy,1)

   if (min_Energy_position == 1) then
     lambda = 0.0D0
     return
   elseif (min_Energy_position == n_points) then
     lambda = step_size * dble(n_points)
     return
   endif

   ! Obtains optimal value parametrizing a quadratic function.
   dE2 = abs(Energy(min_Energy_position) - Energy(min_Energy_position+1))
   dE1 = abs(Energy(min_Energy_position) - Energy(min_Energy_position-1))

   modif_fac = step_size * (dE2 - dE1) / (dE1 + dE2)
   lambda    = step_size * dble(min_Energy_position) - 0.5D0 * modif_fac

end subroutine line_search

subroutine line_search_2d(energy, step_size1, step_size2, lambda1, lambda2)
   ! Performs a linear search for a given Energy array.
   !
   !  If minimum value of E is E(n_points)
   !     return lambda = step_size * n_points
   !  If minimum value of E is E(1)
   !     return lambda = 0.d0
   !  else
   !     returns an expected lambda that minimice E(lambda)
   !     using a parabolic interpolation
   implicit none
   real(kind=8), intent(in)  :: energy(:,:), step_size1, step_size2
   real(kind=8), intent(out) :: lambda1, lambda2

   integer      :: min_energy_position(2)
   real(kind=8) :: dE1, dE2, modif_fac

   if ((size(energy,1)+ size(energy,2))<= 5) then
      write(*,*) "Wrong n_points in lineal search, n_points need to be > 2"
      stop
   end if

   ! Finds minimum value in Energy elements
   min_energy_position = minloc(energy)

   if (min_energy_position(1) == 1) then
     lambda1 = 0.0D0
     call line_search(size(energy,2), energy(:,2), step_size2, lambda2)
     return
   elseif (min_energy_position(1) == size(energy,1)) then
     lambda1 = step_size1 * dble(size(energy,1))
     call line_search(size(energy,2), energy(:,2), step_size2, lambda2)
     return
   endif

   if (min_energy_position(2) == 1) then
      lambda2 = 0.0D0
      call line_search(size(energy,1), energy(:,1), step_size1, lambda1)
      return
    elseif (min_energy_position(2) == size(energy,2)) then
      lambda2 = step_size2 * dble(size(energy,2))
      call line_search(size(energy,1), energy(:,1), step_size1, lambda1)
      return
    endif

   ! Obtains optimal value parametrizing a quadratic function.
   dE1 = abs(Energy(min_energy_position(1)   , min_energy_position(2)) - &
             Energy(min_energy_position(1) +1, min_energy_position(2)))
   dE2 = abs(Energy(min_energy_position(1)   , min_energy_position(2)) - &
             Energy(min_energy_position(1) -1, min_energy_position(2)))
   
   modif_fac = step_size1 * (dE2 - dE1) / (dE1 + dE2)
   lambda1   = step_size1 * dble(min_Energy_position(1)) - 0.5D0 * modif_fac

   dE1 = abs(Energy(min_energy_position(1), min_energy_position(2)   ) - &
             Energy(min_energy_position(1), min_energy_position(2) +1))
   dE2 = abs(Energy(min_energy_position(1), min_energy_position(2)   ) - &
             Energy(min_energy_position(1), min_energy_position(2) -1))

   modif_fac = step_size2 * (dE2 - dE1) / (dE1 + dE2)
   lambda2   = step_size2 * dble(min_Energy_position(2)) - 0.5D0 * modif_fac
end subroutine line_search_2d
