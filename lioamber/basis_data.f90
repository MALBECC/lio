!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  module basis_data
  implicit none
  integer                            :: basis_size
  integer                            :: basis_size_s
  integer                            :: basis_size_p
  integer                            :: basis_size_d
  integer                            :: maximum_contractions
  integer,allocatable,dimension(:)   :: orbital_contractions
  integer,allocatable,dimension(:)   :: parent_atom
  integer,allocatable,dimension(:,:) :: angular_momentum
  real*8,allocatable,dimension(:,:)  :: gauss_expo
  real*8,allocatable,dimension(:,:)  :: gauss_coef
  contains



!------------------------------------------------------------------------------!
  subroutine basis_data_set(ns,np,nd,orba,orbc,ge,gc)
!------------------------------------------------------------------------------!
  implicit none
  integer,intent(in)                :: ns,np,nd
  integer,intent(in),dimension(:)   :: orba
  integer,intent(in),dimension(:)   :: orbc
  real*8,intent(in), dimension(:,:) :: ge
  real*8,intent(in), dimension(:,:) :: gc
  integer                           :: nn


  basis_size_s=ns
  basis_size_p=np
  basis_size_d=nd
!  basis_size_f=nf
  basis_size=ns+np+nd
  maximum_contractions=max(size(ge,2),size(gc,2))

  if (allocated(parent_atom))          deallocate(parent_atom)
  if (allocated(orbital_contractions)) deallocate(orbital_contractions)
  if (allocated(angular_momentum))     deallocate(angular_momentum)
  if (allocated(gauss_expo))           deallocate(gauss_expo)
  if (allocated(gauss_coef))           deallocate(gauss_coef)

  allocate(parent_atom(basis_size))
  allocate(orbital_contractions(basis_size))
  allocate(angular_momentum(3,basis_size))
  allocate(gauss_expo(maximum_contractions,basis_size))
  allocate(gauss_coef(maximum_contractions,basis_size))


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

  return;end subroutine
  end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
