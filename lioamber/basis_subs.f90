!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module basis_subs
   implicit none
   contains
!------------------------------------------------------------------------------!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine basis_data_set(ns,np,nd,orba,orbc,ge,gc)
!------------------------------------------------------------------------------!
   use basis_data, only: basis_size, basis_size_s, basis_size_p, basis_size_d  &
                      &, maximum_contractions, orbital_contractions            &
                      &, parent_atom, angular_momentum, gauss_expo, gauss_coef
   implicit none
   integer,intent(in)                 :: ns,np,nd
   integer,intent(in), dimension(:)   :: orba
   integer,intent(in), dimension(:)   :: orbc
   real*8,intent(in) , dimension(:,:) :: ge
   real*8,intent(in) , dimension(:,:) :: gc
   integer                            :: nn


   basis_size_s=ns
   basis_size_p=np
   basis_size_d=nd
!   basis_size_f=nf
   basis_size=ns+np+nd
   maximum_contractions=max(size(ge,2),size(gc,2))

   if (allocated(parent_atom))          deallocate(parent_atom)
   if (allocated(orbital_contractions)) deallocate(orbital_contractions)
   if (allocated(angular_momentum))     deallocate(angular_momentum)
   if (allocated(gauss_expo))           deallocate(gauss_expo)
   if (allocated(gauss_coef))           deallocate(gauss_coef)

   allocate( parent_atom(basis_size) )
   allocate( orbital_contractions(basis_size) )
   allocate( angular_momentum(3, basis_size) )
   allocate( gauss_expo(maximum_contractions, basis_size) )
   allocate( gauss_coef(maximum_contractions, basis_size) )


   do nn=1,ns
      parent_atom(nn+0)          = orba(nn)
      orbital_contractions(nn+0) = orbc(nn)

      angular_momentum(:,nn) = 0
      gauss_expo(:,nn+0)     = ge(nn,:)
      gauss_coef(:,nn+0)     = gc(nn,:)
   enddo

   do nn=ns+1,ns+np,3
      parent_atom(nn+2)          = orba(nn)
      parent_atom(nn+1)          = orba(nn)
      parent_atom(nn+0)          = orba(nn)
      orbital_contractions(nn+2) = orbc(nn)
      orbital_contractions(nn+1) = orbc(nn)
      orbital_contractions(nn+0) = orbc(nn)

      angular_momentum(:,nn+2) = 0
      angular_momentum(:,nn+1) = 0
      angular_momentum(:,nn+0) = 0
      angular_momentum(1,nn+0) = 1  ! px
      angular_momentum(2,nn+1) = 1  ! py
      angular_momentum(3,nn+2) = 1  ! pz

      gauss_expo(:,nn+2) = ge(nn,:)
      gauss_expo(:,nn+1) = ge(nn,:)
      gauss_expo(:,nn+0) = ge(nn,:)

      gauss_coef(:,nn+2) = gc(nn,:)
      gauss_coef(:,nn+1) = gc(nn,:)
      gauss_coef(:,nn+0) = gc(nn,:)
   enddo

   do nn=ns+np+1,ns+np+nd,6
      parent_atom(nn+5)          = orba(nn)
      parent_atom(nn+4)          = orba(nn)
      parent_atom(nn+3)          = orba(nn)
      parent_atom(nn+2)          = orba(nn)
      parent_atom(nn+1)          = orba(nn)
      parent_atom(nn+0)          = orba(nn)
      orbital_contractions(nn+5) = orbc(nn)
      orbital_contractions(nn+4) = orbc(nn)
      orbital_contractions(nn+3) = orbc(nn)
      orbital_contractions(nn+2) = orbc(nn)
      orbital_contractions(nn+1) = orbc(nn)
      orbital_contractions(nn+0) = orbc(nn)

      angular_momentum(:,nn+5) = 0
      angular_momentum(:,nn+4) = 0
      angular_momentum(:,nn+3) = 0
      angular_momentum(:,nn+2) = 0
      angular_momentum(:,nn+1) = 0
      angular_momentum(:,nn+0) = 0
      angular_momentum(1,nn+0) = 2  ! dxx (x)
      angular_momentum(1,nn+1) = 1  ! dxy (x)
      angular_momentum(2,nn+1) = 1  ! dxy (y)
      angular_momentum(2,nn+2) = 2  ! dyy (y)
      angular_momentum(1,nn+3) = 1  ! dxz (x)
      angular_momentum(3,nn+3) = 1  ! dxz (z)
      angular_momentum(2,nn+4) = 1  ! dyz (y)
      angular_momentum(3,nn+4) = 1  ! dyz (z)
      angular_momentum(3,nn+5) = 2  ! dzz (z)

      gauss_expo(:,nn+5) = ge(nn,:)
      gauss_expo(:,nn+4) = ge(nn,:)
      gauss_expo(:,nn+3) = ge(nn,:)
      gauss_expo(:,nn+2) = ge(nn,:)
      gauss_expo(:,nn+1) = ge(nn,:)
      gauss_expo(:,nn+0) = ge(nn,:)

      gauss_coef(:,nn+5) = gc(nn,:)/SQRT(3.0)
      gauss_coef(:,nn+4) = gc(nn,:)
      gauss_coef(:,nn+3) = gc(nn,:)
      gauss_coef(:,nn+2) = gc(nn,:)/SQRT(3.0)
      gauss_coef(:,nn+1) = gc(nn,:)
      gauss_coef(:,nn+0) = gc(nn,:)/SQRT(3.0)
   enddo

end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine basis_data_norm( Isize, Icont, gcoefs )
!------------------------------------------------------------------------------!
   use maskrmm,     only: rmmget_fock
   use faint_cpu,   only: int1
   use basis_data,  only: gauss_coef
   use garcha_mod,  only: RMM,Nuc,a,c,d,r,Iz,NORM,natom,M,Md,ncont,nshell,ntatom

   implicit none
   integer, intent(in)              :: Isize
   integer, intent(in)              :: Icont
   double precision , intent(inout) :: gcoefs(Isize, Icont)

   integer                          :: kk, MM, M5, M11
   double precision                 :: scratch_energy
   double precision, allocatable    :: Smat(:,:)

!   call g2g_timer_start('RMMcalc1')
   allocate( Smat(isize, isize) )
   MM = M*(M+1)/2 ; M5 = 1 + MM*2; M11 = M5+MM+Md*(Md+1)
   call int1(scratch_energy, RMM(M5:M5+MM), RMM(M11:M11+MM), Smat, d, r, Iz, &
             natom, ntatom)
   do kk = 1, isize
      gcoefs(kk,:) = gcoefs(kk,:) / sqrt( Smat(kk,kk) )
      gauss_coef(:,kk) = gauss_coef(:,kk) / sqrt( Smat(kk,kk) )
   enddo
   deallocate( Smat )
!   call g2g_timer_stop('RMMcalc1')

end subroutine basis_data_norm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
