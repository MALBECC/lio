subroutine reduced_space(Cin,Ein,Cout,Eout,map_occ,map_vir,NCO,M,NCOlr,Mlr,Nvirt,Ndim)
! This routine perform the Reduced Subspace in Linear Response according
! Besley' paper. DOI: 10.1016/j.cplett.2004.04.004
! and Mulliken's paper. DOI: 10.1063/1.1740588

! The reduced.dat file needed for this approximation must
! contain 2 lines
! 1) how many atoms
! 2) the id of those atoms
   implicit none

   integer, intent(in) :: NCO, M
   integer, intent(out) :: NCOlr, Mlr, Nvirt, Ndim
   LIODBLE, intent(in) :: Cin(:,:), Ein(:)
   integer, allocatable, intent(out) :: map_occ(:), map_vir(:)
   LIODBLE, allocatable, intent(out) :: Cout(:,:), Eout(:)

   logical :: res
   integer :: ii, red_natom, ind
   integer, dimension(:)  , allocatable :: map_temp, red_list
   LIODBLE, dimension(:,:), allocatable :: Mmo

   print*, "*Using Reduced atoms MOs Sub-space"

   ! Read atoms of subspace
   res = .false.
   inquire(file="reduced.dat",exist=res)
   if ( res ) then
      open(unit=456,file="reduced.dat")
      read(456,*) red_natom
      allocate(red_list(red_natom))
      read(456,*) (red_list(ii), ii=1,red_natom)
      close(456)
   else
      print*, "The file reduced.dat does not exist"
      stop
   endif

   print*, "The atoms considered are"
   do ii=1,red_natom
      write(*,"(1X,I3)",ADVANCE="NO") red_list(ii)
   enddo
   write(*,*)
   
!  We calculate the Mulliken matrix in AO,MO
   allocate(Mmo(M,NCO))
   call get_Mmo(Cin,Mmo,M,NCO)

!  Obtain Reduced Subspace: Occupied
   allocate(map_temp(NCO)); map_temp = 0
   call get_occupied(Mmo,red_list,map_temp,NCO,red_natom,M,NCOlr)
   if(allocated(map_occ)) deallocate(map_occ)
     allocate(map_occ(NCOlr)) ! map_occ(ind small) -> ind big
   map_occ(1:NCOlr) = map_temp(1:NCOlr)
   deallocate(map_temp,Mmo)

!  Obtain Reduced Subspace: Virtual
   allocate(map_temp(M-NCO)); map_temp = 0
   call get_virtual(Cin,red_list,M,NCO,red_natom,Nvirt,map_temp)
   if(allocated(map_vir)) deallocate(map_vir)
     allocate(map_vir(Nvirt))
   map_vir(1:Nvirt) = map_temp(1:Nvirt)
   deallocate(map_temp)

!  Get new dimensions
   Mlr  = NCOlr + Nvirt
   Ndim = NCOlr * Nvirt
   
   if(allocated(Cout)) deallocate(Cout)
   if(allocated(Eout)) deallocate(Eout)
     allocate(Cout(M,Mlr),Eout(Mlr))

!  Get MOs and Energies truncated
   do ii=1,NCOlr
      ind = map_occ(ii)
      Cout(:,ii) = Cin(:,ind)
      Eout(ii) = Ein(ind)
   enddo
   do ii=1,Nvirt
      ind = NCO + map_vir(ii)
      Cout(:,NCOlr+ii) = Cin(:,ind)
      Eout(NCOlr+ii) = Ein(ind)
   enddo
   print*, " "

   deallocate(red_list)
end subroutine reduced_space

subroutine get_Mmo(Cin,Mmo,M,NCO)
use garcha_mod, only: Iz, d, r, natom, ntatom
use basis_data, only: MM
use faint_cpu , only: int1
   implicit none

   integer, intent(in) :: M, NCO
   LIODBLE, intent(in) :: Cin(:,:)
   LIODBLE, intent(out):: Mmo(M,NCO)

   integer :: ii, jj, kk
   LIODBLE :: En, temp
   LIODBLE, dimension(:,:), allocatable :: Sao
   LIODBLE, dimension(:)  , allocatable :: F1e, Hmat

   allocate(Sao(M,M),F1e(MM),Hmat(MM)) 

   call int1(En, F1e, Hmat, Sao, d, r, Iz, natom, ntatom)
   call spunpack('L', M, F1e, Sao)
   do ii=1,NCO
      do jj=1,M
         temp = 2.0d0 * Cin(jj,ii)*Cin(jj,ii)
         do kk=jj+1,M
            temp = temp + 4.0d0 * Cin(jj,ii)*Cin(kk,ii)*Sao(jj,kk)
         enddo
         Mmo(jj,ii) = temp
      enddo
   enddo

   deallocate(Sao,F1e,Hmat)

end subroutine get_Mmo

subroutine get_occupied(Mmo,red_list,map,NCO,red_natom,M,NCOlr)
use basis_data  , only: Nuc
use excited_data, only: thres_occ
   implicit none

   integer, intent(in) :: NCO, red_natom, M
   integer, intent(in) :: red_list(red_natom)
   LIODBLE, intent(in) :: Mmo(M,NCO)
   integer, intent(out):: NCOlr, map(NCO)

   integer :: ii, jj, kk, cant
   LIODBLE :: temp

   write(*,"(1X,A,F8.6)") "occ threshold ", thres_occ
   cant = 0
   do ii=1,NCO
      temp = 0.0d0
      do jj=1,M
         do kk=1,red_natom
            if ( Nuc(jj) == red_list(kk) ) then
                temp = temp + Mmo(jj,ii)
            endif
         enddo
      enddo
      if ( temp >= thres_occ ) then
         cant = cant + 1
         map(cant) = ii
      endif
   enddo
   NCOlr = cant

end subroutine get_occupied

subroutine get_virtual(Cin,red_list,M,NCO,red_natom,Nvirt,map)
use basis_data  , only: Nuc
use excited_data, only: thres_vir
   implicit none

   integer, intent(in) :: M, NCO, red_natom
   integer, intent(in) :: red_list(red_natom)
   LIODBLE, intent(in) :: Cin(:,:)
   integer, intent(out):: Nvirt, map(M-NCO)

   integer :: ii, jj, kk, ind, cant
   LIODBLE :: temp
   LIODBLE, dimension(:,:), allocatable :: Cvir

   write(*,"(1X,A,F8.6)") "vir threshold ", thres_vir

   ! Normalization of Virtual molecular orbitals
   allocate(Cvir(M,M-NCO))
   do ii=1,M-NCO
      temp = 0.0d0
      ind = NCO + ii
      do jj=1,M
         temp = temp + Cin(jj,ind)*Cin(jj,ind)
      enddo
      Cvir(:,ii) = Cin(:,ind) / dsqrt(temp)
   enddo

   cant = 0
   do ii=1,M-NCO
      temp = 0.0d0
      do jj=1,M
         do kk=1,red_natom
            if ( Nuc(jj) == red_list(kk) ) then
                temp = temp + Cvir(jj,ii)*Cvir(jj,ii)
            endif
         enddo
      enddo

      if ( temp >= thres_vir ) then
         cant = cant + 1
         map(cant) = ii
      endif

   enddo
   Nvirt = cant

   deallocate(Cvir)
end subroutine get_virtual
