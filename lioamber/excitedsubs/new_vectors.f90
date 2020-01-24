subroutine new_vectors(R,Esub,Eall,T,Eold,Ndim,Sdim,&
                       Nstat,Mlr,NCO,New,conv)
use excited_data, only: tolv, tole
   implicit none

   integer, intent(in) :: Ndim, Sdim, Nstat, Mlr, NCO
   double precision, intent(in) :: R(Ndim,Nstat), Esub(Nstat), Eall(Mlr)
   double precision, intent(inout) :: T(Ndim,Sdim+Nstat), Eold(Nstat)
   integer, intent(out) :: New
   logical, intent(out) :: conv

   integer :: i, iv, occ, virt, Nvirt, ind, NCOc, append 
   double precision :: temp, ERROR, MAX_ERROR, MAX_ENE, diffE, norm
   double precision, dimension(:,:), allocatable :: Qvec
   integer, dimension(:), allocatable :: valid_id

   allocate(Qvec(Ndim,Nstat))
   allocate(valid_id(Nstat))

   MAX_ERROR = 0.0D0
   MAX_ENE = 0.0D0
   Nvirt = Mlr - NCO
   NCOc = NCO + 1

   do iv=1,Nstat
     diffE = abs(Eold(iv) - Esub(iv))
     if(diffE > MAX_ENE) MAX_ENE = diffE

     call norma(R(:,iv),Ndim,ERROR)
     if(ERROR > MAX_ERROR) MAX_ERROR = ERROR

     if(ERROR > tolv .or. diffE > tole) then
        New = New + 1
        valid_id(New) = iv
        do i=1,Ndim
          ind = i - 1
          occ = NCO - (ind/Nvirt)
          virt = mod(ind,Nvirt) + NCOc
          temp = Eall(virt) - Eall(occ)
          temp = 1.0D0 / (Esub(iv) - temp)
          Qvec(i,iv) = temp * R(i,iv)
        enddo
     else
        write(*,"(1X,A,I2,1X,A)") "Vector:",iv,"Converged"
     endif
   enddo

   if(Sdim + New >= Ndim) then
      conv = .true.
      return
   endif

   write(*,8070) MAX_ERROR, tolv, MAX_ENE, tole
   if(MAX_ERROR < tolv .and. MAX_ENE < tole) then
     conv = .true.
   else ! Append New vectors
     conv = .false.
     append = 0

     ! Orthonormalization
     do iv=1,New
        do i=1,Sdim+append
           call prod(T(1:Ndim,i),Qvec(1:Ndim,valid_id(iv)),norm,Ndim)
           Qvec(1:Ndim,valid_id(iv))=Qvec(1:Ndim,valid_id(iv)) - norm * T(1:Ndim,i)
        enddo
        call norma(Qvec(1:Ndim,valid_id(iv)),Ndim,ERROR)
        ERROR = 1.0d0 / dsqrt(ERROR)
        T(1:Ndim,Sdim+iv) = ERROR * Qvec(1:Ndim,valid_id(iv))
        append = append + 1
     enddo
   endif

   Eold = Esub
   deallocate(Qvec,valid_id)

8070   FORMAT(1X,"VectorsError (crit) = ", F15.7," (",ES9.2,")" &
              " - EnergyError (crit) = ", F15.7," (",ES9.2,")" )
end subroutine new_vectors

subroutine prod(A,B,temp,N)
   implicit none

   integer, intent(in) :: N
   double precision, intent(in) :: A(N), B(N)
   double precision, intent(out) :: temp

   integer :: ii

   temp = 0.0D0

   do ii=1,N
      temp = temp + A(ii) * B(ii)
   enddo
end subroutine prod
