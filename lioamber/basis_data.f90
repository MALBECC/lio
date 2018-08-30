!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% BASIS.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains two modules: basis_data, containing basis set function    !
! data, and basis_subs, containing the following subroutines:                  !
!                                                                              !
! External access:                                                             !
! · basis_init()       : Initializes basis data.                               !
!    --> Requires basis and fitting sets names (or basis file name), number of !
!        atoms and nuclear atomic charges.                                     !
! · basis_setup_ehren(): Initialises data for ehrenfest runs                   !
!    --> No inputs required.                                                   !
! · basis_deinit()     : Deallocates basis data.                               !
!    --> No inputs required.                                                   !
!                                                                              !
! Intended for basis_subs internal-only access:                                !
! · basis_read_internal(): Reads basis functions from /dat directory.          !
! · basis_read_external(): Reads basis functions form a custom basis file.     !
! · basis_set_size()     : Prereads data and sets array sizes                  !
! · check_basis()        : Checks if all atoms have basis functions  assigned. !
!                                                                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

! TODO: INVERT C and A matrices, since loops are faster that way. In Ehrenfest
!       these are already inverted.

module basis_data
   implicit none
   ! Namelist inputfiles
   ! basis_set  : The name of the basis set or basis_set file.
   ! fitting_set: The name of the fitting set.
   ! int_basis  : Use LIO internal basis (located in /dat)
   ! rMax       : Maximum exponent for double-precision integrals.
   ! rMaxs      : Maximum exponent for single-precision integrals.
   ! norm       : Normalize integrals (deprecated).
   character(len=100) :: basis_set   = "DZVP"
   character(len=100) :: fitting_set = "DZVP Coulomb Fitting"
   logical            :: int_basis   = .false.
   double precision   :: rMax        = 16.0D0
   double precision   :: rMaxs       =  5.0D0
   logical            :: norm        = .true.

   ! Single variables
   ! M     : number of basis functions.
   ! Md    : number of fitting functions.
   ! MM    : M*(M+1), used for vectorized matrices.
   ! MMd   : Md*(Md+1), used for vectorized matrices.
   ! kknums: number of single-precision two-center integrals.
   ! kknumd: number of double-precision two-center integrals.
   integer          :: M      = 0
   integer          :: Md     = 0
   integer          :: MM     = 0
   integer          :: MMd    = 0
   integer          :: kknums = 0
   integer          :: kknumd = 0

   ! Preset arrays
   ! nShell : Number of basis functions for each shell (s,p,d,f)
   ! nShelld: Number of auxiliary basis functions for each shell (s,p,d,f)
   ! max_c_per_atom: Maximum number of contractions for a single function.
   ! max_f_per_atom: Maximum number of functions for a single atom.
   integer :: nShell(0:3)
   integer :: nShelld(0:3)
   integer :: max_c_per_atom
   integer :: max_f_per_atom

   ! Allocatable arrays
   ! Nuc(i)     : Index of the atom containing function i.
   ! Nucd(i)    : Index of the atom containing fitting function i.
   ! nCont(i)   : Number of contractions for function i.
   ! nContd(i)  : Number of contractions for fitting function i.
   ! ang_mom(i) : Angular momentum of function i. Or an angry mother. Your choice.
   ! ang_momd(i): Angular momentum of auxiliary function i.
   ! a(i,j)     : Exponent for function i, contraction j.
   ! c(i,j)     : Coefficient for function i, contraction j.
   ! ad(i,j)    : Exponent for auxiliary function i, contraction j.
   ! cd(i,j)    : Coefficient for auxiliary function i, contraction j.
   ! atmin(i)   : The minimum exponent found for atom i.
   ! indexii(i) : Function index after reordering by s, p, d.
   ! indexiid(i): Auxiliary function index after reordering by s, p, d.
   ! kkInd(:)   : Index for double-precision two-center integrals.
   ! kkInds(:)  : Index for single-precision two-center integrals.
   ! cool(:)    : Storage for two-center integrals in  double precision.
   ! cools(:)   : Storage for two-center integrals in  single precision.
   integer         , allocatable :: Nuc(:)
   integer         , allocatable :: Nucd(:)
   integer         , allocatable :: nCont(:)
   integer         , allocatable :: nContd(:)
   integer         , allocatable :: ang_mom(:)
   integer         , allocatable :: ang_momd(:)
   integer         , allocatable :: indexii(:)
   integer         , allocatable :: indexiid(:)
   integer         , allocatable :: kkInd(:)
   integer         , allocatable :: kkInds(:)
   integer         , allocatable :: natomc(:)
   integer         , allocatable :: nns(:)
   integer         , allocatable :: nnp(:)
   integer         , allocatable :: nnd(:)
   integer         , allocatable :: nnps(:)
   integer         , allocatable :: nnpd(:)
   integer         , allocatable :: nnpp(:)
   integer         , allocatable :: jatc(:,:)
   double precision, allocatable :: a(:,:)
   double precision, allocatable :: c(:,:)
   double precision, allocatable :: ad(:,:)
   double precision, allocatable :: cd(:,:)
   double precision, allocatable :: af(:)
   double precision, allocatable :: atmin(:)
   double precision, allocatable :: cool(:)
   real            , allocatable :: cools(:)

   ! Temporary for EHRENFEST
   double precision, allocatable :: a_ehren(:,:)
   double precision, allocatable :: c_ehren(:,:)
   integer         , allocatable :: ang_mom_ehren(:,:)

   ! Temporary for Linear Response
   double precision, allocatable :: c_raw(:,:)

   ! GLOBAL PARAMETERS
   ! Degeneracy for each angular momentum
   integer         , parameter :: ANG_DEG(0:3) = (/1, 3, 6, 10/)
   integer         , parameter :: MAX_CONTRACT = 50
contains
end module basis_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module basis_subs
   implicit none
   integer, parameter :: TMP_OPEN_UID = 9999
contains

subroutine basis_init(basis_name, fitting_name, n_atoms, atom_Z, out_stat)
   use basis_data, only: M, Md, int_basis, Nuc, Nucd, nCont, nContd, a, c, ad, &
                         cd, atmin, nns, nnp, nnd, nshell, nshelld, norm, af,  &
                         indexii, indexiid, natomc, jatc, nnps, nnpp, nnpd,    &
                         ang_mom, ang_momd, max_f_per_atom, max_c_per_atom, MM,&
                         MMd, c_raw
   implicit none
   ! Inputs:
   !   n_atoms        : the number of atoms in the QM system.
   !   atom_Z(n_atoms): nuclear atomic charge of said atoms.
   !   basis_name     : the name of the basis set or basis file.
   !   fitting_name   : the name of the fitting set (not used when int_basis is
   !                    true).
   ! Outputs:
   !   out_stat       : Returns an iostat-like value for error handling.

   integer         , intent(in)    :: n_atoms, atom_Z(n_atoms)
   integer         , intent(out)   :: out_stat
   character(len=*), intent(inout) :: basis_name, fitting_name

   integer :: icount, iostat
   integer, allocatable :: atom_count(:)
   logical, allocatable :: atom_basis_chk(:), atom_fitting_chk(:)

   out_stat = 0
   allocate(atom_count(0:120), atom_basis_chk(0:120), atom_fitting_chk(0:120))

   ! Precalculates the amount of atoms of a certain type, in order to
   ! correctly asses the sizes of M and Md
   atom_count = 0
   do icount = 1, n_atoms
      atom_count(atom_Z(icount)) = atom_count(atom_Z(icount)) + 1
   enddo

   call basis_set_size(M, Md, max_f_per_atom, max_c_per_atom, atom_basis_chk,  &
                       atom_fitting_chk, basis_name, fitting_name, atom_count, &
                       atom_Z, n_atoms, int_basis, iostat)
   if (iostat .gt. 0) then
      out_stat = 1
      return
   endif

   call check_basis(atom_basis_chk, atom_fitting_chk, n_atoms, atom_Z, iostat)
   if (iostat .gt. 0) then
      out_stat = 2
      return
   endif

   allocate(c(M, max_c_per_atom), a(M, max_c_per_atom), cd(Md, max_c_per_atom),&
            ad(Md, max_c_per_atom), nCont(M), nContd(Md), ang_mom(M),          &
            ang_momd(Md), Nuc(M), Nucd(Md), indexii(M), indexiid(Md), af(Md),  &
            c_raw(M, max_c_per_atom))
   allocate(atmin(n_atoms), nns(n_atoms), nnp(n_atoms), nnd(n_atoms),     &
            nnps(n_atoms), nnpp(n_atoms), nnpd(n_atoms), natomc(n_atoms), &
            jatc(n_atoms,n_atoms))

   ! Initializes everything to 0.
   c  = 0.0D0 ; a  = 0.0D0; nCont  = 0; ang_mom  = 0; nuc  = 0; af = 0.0D0
   cd = 0.0D0 ; ad = 0.0D0; nContd = 0; ang_momd = 0; nucd = 0

   nns  = 0 ; nnp  = 0  ; nnd  = 0 ; natomc = 0
   nnps = 0 ; nnpp = 0  ; nnpd = 0 ; jatc   = 0


   if (int_basis) then
      call read_basis_internal(basis_name, fitting_name, M, Md, n_atoms, norm, &
                               max_f_per_atom, max_c_per_atom, atom_Z, c_raw,  &
                               c, a, cd, ad, nCont, nContd, ang_mom, ang_momd, &
                               Nuc, Nucd, atmin, nns, nnp, nnd, nShell,        &
                               nShelld, iostat)
   else
      call read_basis_external(basis_name, M,Md, n_atoms, norm, max_f_per_atom,&
                               max_c_per_atom, atom_Z, c_raw, c, a, cd, ad,    &
                               nCont, nContd, ang_mom, ang_momd, Nuc, Nucd,    &
                               atmin, nns, nnp, nnd, nShell, nShelld, iostat)
   endif
   if (iostat .gt. 0) then
      out_stat = 3
      return
   endif

   ! Reorders basis: first all s, then all p, then all d.
   call reorder_basis(a,  c,  Nuc,  nCont,  indexii,  M,  max_c_per_atom, &
                      ang_mom,  nShell, c_raw)
   call reorder_basis(ad, cd, Nucd, nContd, indexiid, Md, max_c_per_atom, &
                      ang_momd, nShelld)

   ! Sets MM and MMd
   MM = M *(M +1) / 2 ; MMd = Md *(Md +1) / 2 

   deallocate(atom_count, atom_basis_chk, atom_fitting_chk, ang_mom, ang_momd)
end subroutine basis_init

subroutine basis_deinit()
   use basis_data, only: Nuc, Nucd, nCont, nContd, a, c, ad, cd, atmin, nns, &
                         nnp, nnd, af, indexii, indexiid, natomc, jatc, nnps,&
                         nnpp, nnpd, c_raw

   implicit none

   ! M or Md sized.
   if (allocated(c))        deallocate(c)
   if (allocated(a))        deallocate(a)
   if (allocated(cd))       deallocate(cd)
   if (allocated(ad))       deallocate(ad)
   if (allocated(af))       deallocate(af)
   if (allocated(ncont))    deallocate(ncont)
   if (allocated(ncontd))   deallocate(ncontd)
   if (allocated(nuc))      deallocate(nuc)
   if (allocated(nucd))     deallocate(nucd)
   if (allocated(indexii))  deallocate(indexii)
   if (allocated(indexiid)) deallocate(indexiid)

   ! natom sized.
   if (allocated(atmin))  deallocate(atmin)
   if (allocated(jatc))   deallocate(jatc)
   if (allocated(natomc)) deallocate(natomc)
   if (allocated(nns))    deallocate(nns)
   if (allocated(nnp))    deallocate(nnp)
   if (allocated(nnd))    deallocate(nnd)
   if (allocated(nnps))   deallocate(nnps)
   if (allocated(nnpp))   deallocate(nnpp)
   if (allocated(nnpd))   deallocate(nnpd)


   ! Ehrenfest and LR-TDDFT
   if (allocated(c_raw)) deallocate(c_raw)

end subroutine basis_deinit

subroutine basis_setup_ehren()
   use basis_data, only: a_ehren, c_ehren, ang_mom_ehren, max_c_per_atom, &
                         a, c, nShell, M

   implicit none
   integer :: icount

   ! Transposes coefficient and exponent matrices.
   allocate (c_ehren(max_c_per_atom,M), a_ehren(max_c_per_atom, M))
   do icount = 1, M
      a_ehren(:,icount) = a(icount,:)
      c_ehren(:,icount) = c(icount,:)
   enddo

   ! Calculates angular momenta and directions
   allocate(ang_mom_ehren(3,M))
   ang_mom_ehren = 0
   do icount = nShell(0)+1, nShell(0)+nShell(1), 3
      ang_mom_ehren(1, icount  ) = 1
      ang_mom_ehren(2, icount+1) = 1
      ang_mom_ehren(3, icount+2) = 1
   enddo
   do icount = nShell(0)+nShell(1)+1, nShell(0)+nShell(1)+nShell(2), 6
      ang_mom_ehren(1,icount+0) = 2  ! dxx (x)
      ang_mom_ehren(1,icount+1) = 1  ! dxy (x)
      ang_mom_ehren(2,icount+1) = 1  ! dxy (y)
      ang_mom_ehren(2,icount+2) = 2  ! dyy (y)
      ang_mom_ehren(1,icount+3) = 1  ! dxz (x)
      ang_mom_ehren(3,icount+3) = 1  ! dxz (z)
      ang_mom_ehren(2,icount+4) = 1  ! dyz (y)
      ang_mom_ehren(3,icount+4) = 1  ! dyz (z)
      ang_mom_ehren(3,icount+5) = 2  ! dzz (z)

      c_ehren(:,icount  ) = c_ehren(:,icount  ) / dsqrt(3.0D0)
      c_ehren(:,icount+2) = c_ehren(:,icount+2) / dsqrt(3.0D0)
      c_ehren(:,icount+5) = c_ehren(:,icount+5) / dsqrt(3.0D0)
   enddo
end subroutine basis_setup_ehren

subroutine basis_set_size(basis_size, aux_size, max_f_per_atom, max_c_per_atom,&
                          atom_bas_done, atom_fit_done, basis_file,            &
                          fitting_file, atom_count, atom_Z, n_atoms,           &
                          use_internal, iostatus)
   use basis_data, only: ANG_DEG, MAX_CONTRACT
   implicit none
   integer         , intent(in)     :: n_atoms, atom_Z(n_atoms), &
                                       atom_count(0:120)
   logical         , intent(in)     :: use_internal
   integer         , intent(out)    :: basis_size, aux_size, max_f_per_atom, &
                                       max_c_per_atom, iostatus
   logical         , intent(out)    :: atom_bas_done(0:120), &
                                       atom_fit_done(0:120)
   character(len=*), intent(inout)  :: basis_file, fitting_file

   logical              :: file_exists
   integer              :: file_uid = TMP_OPEN_UID
   integer              :: file_iostat, icount, iatom, nraw, ncon
   character(len=20)    :: start_str
   character(len=100)   :: lio_dir, line_read
   integer, allocatable :: l_of_func(:)

   allocate(l_of_func(MAX_CONTRACT))
   l_of_func  = 0
   iostatus   = 0
   basis_size = 0
   aux_size   = 0

   if (.not. use_internal) then
      ! Read basis from external input basis file.
      inquire(file = basis_file, exist = file_exists)
      if (.not. file_exists) then
         write(*,'(A)') "  Error: Basis set file ", trim(basis_file), &
                        " not found."
         iostatus = 1
         return
      endif
      open(unit = file_uid, file= basis_file, iostat = file_iostat)
      read(file_uid, '(A8)') start_str

      do while (start_str .ne. "endbasis")
         read(file_uid,*) iatom, nraw, ncon
         if (max_f_per_atom .lt. nraw) max_f_per_atom = nraw
         if (max_c_per_atom .lt. ncon) max_c_per_atom = ncon
         if (ncon .gt. MAX_CONTRACT) then
            write(*,'(A)') "  Error: Atom has more contractions than the &
                           &maximum allowed (", MAX_CONTRACT,")."
            iostatus = 2
            return
         endif

         ! Only needs to read angular momenta for each contraction, the rest is
         ! skipped.
         read(file_uid,*)
         read(file_uid,*) (l_of_func(icount), icount=1, ncon)
         if (any(atom_Z == iatom)) then
            do icount = 1, ncon
               basis_size = basis_size &
                          + ANG_DEG(l_of_func(icount)) * atom_count(iatom)
            enddo
         endif
         do icount = 1, nraw
            read(file_uid,*)
         enddo
         atom_bas_done(iatom) = .true.

         ! Reads auxiliary basis set for an atom doing the same as before.
         read(file_uid,*) iatom, nraw, ncon
         if (max_f_per_atom .lt. nraw) max_f_per_atom = nraw
         if (max_c_per_atom .lt. ncon) max_c_per_atom = ncon
         if (ncon .gt. MAX_CONTRACT) then
            write(*,'(A)') "  Error: Atom has more contractions than the &
                           &maximum allowed (", MAX_CONTRACT,")."
            iostatus = 2
            return
         endif

         read(file_uid,*)
         read(file_uid,*) (l_of_func(icount), icount=1, ncon)
         if (any(atom_Z == iatom)) then
            do icount = 1, ncon
               aux_size = aux_size &
                          + ANG_DEG(l_of_func(icount)) * atom_count(iatom)
            enddo
         endif
         do icount = 1, nraw
            read(file_uid,*)
         enddo
         atom_fit_done(iatom) = .true.

         read(file_uid, '(A8)') start_str
      enddo
      close(file_uid)
   else
      ! Reads from LIO internal basis set files.
      call getenv("LIOHOME", lio_dir)
      if (lio_dir == "") then
         write(*,'(A)') "  Error: LIOHOME not set. Cannot use internal basis", &
                        " files. Either set LIOHOME to your LIO installation", &
                        " directory or use an external basis set file."
         iostatus = 2
         return
      endif
      basis_file   = trim(lio_dir) // "/dat/basis/" // basis_file
      fitting_file = trim(lio_dir) // "/dat/basis/fitting/" // fitting_file

      ! Checks file existence.
      inquire(file = basis_file, exist = file_exists)
      if (.not. file_exists) then
         write(*,'(A)') "  Error: Basis set file ", trim(basis_file), &
                        " not found."
         iostatus = 1
         return
      endif
      inquire(file = fitting_file, exist = file_exists)
      if (.not. file_exists) then
         write(*,'(A)') "  Error: Fitting set file ", trim(fitting_file), &
                        " not found."
         iostatus = 1
         return
      endif

      ! Reads basis set data.
      ! Skips empty lines and those starting with # (comments)
      open(unit = file_uid, file= basis_file, iostat = file_iostat)
      line_read = ""
      do while (line_read == "")
         read(file_uid,*) line_read
         read(line_read, '(A1)') start_str
         if (start_str == "#") line_read = ""
      enddo
      read(line_read, '(A8)') start_str

      do while (start_str .ne. "endbasis")
         read(file_uid,*) iatom, nraw, ncon
         if (max_f_per_atom .lt. nraw) max_f_per_atom = nraw
         if (max_c_per_atom .lt. ncon) max_c_per_atom = ncon
         if (ncon .gt. MAX_CONTRACT) then
            write(*,'(A)') "  Error: Atom has more contractions than the &
                           &maximum allowed (", MAX_CONTRACT,")."
            iostatus = 2
            return
         endif

         ! Only needs to read angular momenta for each contraction, the rest is
         ! skipped.
         read(file_uid,*)
         read(file_uid,*) (l_of_func(icount), icount=1, ncon)
         if (any(atom_Z == iatom)) then
            do icount = 1, ncon
               basis_size = basis_size &
                          + ANG_DEG(l_of_func(icount)) * atom_count(iatom)
            enddo
         endif
         do icount = 1, nraw
            read(file_uid,*)
         enddo
         atom_bas_done(iatom) = .true.

         ! Skips empty lines or those starting with #.
         line_read = ""
         do while (line_read == "")
            read(file_uid,*) line_read
            read(line_read, '(A1)') start_str
            if (start_str == "#") line_read = ""
         enddo
         read(line_read, '(A8)') start_str
      enddo
      close(file_uid)

      ! Reads fitting set data.
      ! Skips empty lines or those starting with #.
      open(unit = file_uid, file= fitting_file, iostat = file_iostat)
      line_read = ""
      do while (line_read == "")
         read(file_uid,*) line_read
         read(line_read, '(A1)') start_str
         if (trim(start_str) == "#") line_read = ""
      enddo
      read(line_read, '(A8)') start_str

      do while (start_str .ne. "endbasis")
         read(file_uid,*) iatom, nraw, ncon
         if (max_f_per_atom .lt. nraw) max_f_per_atom = nraw
         if (max_c_per_atom .lt. ncon) max_c_per_atom = ncon
         if (ncon .gt. MAX_CONTRACT) then
            write(*,'(A)') "  Error: Atom has more contractions than the &
                           &maximum allowed (", MAX_CONTRACT,")."
            iostatus = 2
            return
         endif

         ! Only needs to read angular momenta for each contraction, the rest is
         ! skipped.
         read(file_uid,*)
         read(file_uid,*) (l_of_func(icount), icount=1, ncon)
         if (any(atom_Z == iatom)) then
            do icount = 1, ncon
               aux_size = aux_size &
                          + ANG_DEG(l_of_func(icount)) * atom_count(iatom)
            enddo
         endif
         do icount = 1, nraw
            read(file_uid,*)
         enddo
         atom_fit_done(iatom) = .true.

         ! Skips empty lines or those starting with #.
         line_read = ""
         do while (line_read == "")
            read(file_uid,*) line_read
            read(line_read, '(A1)') start_str
            if (trim(start_str) == "#") line_read = ""
         enddo
         read(line_read, '(A8)') start_str
      enddo
      close(file_uid)
   endif
   deallocate(l_of_func)
end subroutine basis_set_size

subroutine check_basis(atom_bas_done, atom_fit_done, n_atoms, atom_Z, iostatus)
   use fileio, only: atom_name

   implicit none
   integer, intent(in)  :: n_atoms, atom_Z(n_atoms)
   logical, intent(in)  :: atom_bas_done(0:120), atom_fit_done(0:120)
   integer, intent(out) :: iostatus

   integer          :: iatom
   character(len=3) :: iatom_name

   iostatus = 0
   do iatom = 1, n_atoms
      if ( .not. atom_bas_done(atom_Z(iatom)) ) then
         call atom_name(iatom, iatom_name)
         write(*,'(A)') "  Error: Basis set not found for ", &
                        trim(iatom_name), "."
         iostatus = 1
         return
      endif
      if ( .not. atom_fit_done(atom_Z(iatom)) ) then
         call atom_name(iatom, iatom_name)
         write(*,'(A)') "  Error: Fitting set not found for ", &
                        trim(iatom_name), "."
         iostatus = 1
         return
      endif
   enddo
end subroutine check_basis

subroutine read_basis_external(basis_file, n_funcs, n_fits, n_atoms, normalize,&
                               max_fun_per_atom, max_con_per_atom, atom_Z,     &
                               craw, coef, expo, coefd, expod, n_cont, n_contd,&
                               ang_mom_f, ang_mom_fd, atm_of_func,             &
                               atm_of_funcd, min_atm_exp, nns, nnp, nnd,       &
                               nShell, nShelld, iostatus)
   use basis_data   , only: ANG_DEG
   use constants_mod, only: PI32
   implicit none
   integer         , intent(in)  :: max_con_per_atom, max_fun_per_atom, &
                                    n_atoms, atom_Z(n_atoms), n_funcs, n_fits
   logical         , intent(in)  :: normalize
   character(len=*), intent(in)  :: basis_file
   integer         , intent(out) :: atm_of_func(n_funcs), atm_of_funcd(n_fits),&
                                    nns(n_atoms), nnp(n_atoms), nnd(n_atoms),  &
                                    iostatus, n_cont(n_funcs), n_contd(n_fits),&
                                    ang_mom_f(n_funcs),  ang_mom_fd(n_fits),   &
                                    nShell(0:3), nShelld(0:3)
   double precision, intent(out) :: coef(n_funcs,max_con_per_atom), &
                                    craw(n_funcs,max_con_per_atom), &
                                    expo(n_funcs,max_con_per_atom), &
                                    coefd(n_fits,max_con_per_atom), &
                                    expod(n_fits,max_con_per_atom), &
                                    min_atm_exp(n_atoms)

   integer            :: file_iostat, file_uid = TMP_OPEN_UID
   integer            :: iatom, nraw, ncon, atom, icont, icount, l2, index,    &
                         n_orig, n_aux
   double precision   :: min_exp
   character(len=20)  :: start_str

   logical         , allocatable :: basis_done(:), fitting_done(:)
   integer         , allocatable :: n_cont_func(:), ang_mom(:)
   double precision, allocatable :: expo_temp(:), coef_temp(:)

   allocate(n_cont_func(max_con_per_atom)   , ang_mom(max_con_per_atom)   , &
            expo_temp(max_fun_per_atom)     , coef_temp(max_fun_per_atom) , &
            basis_done(n_atoms)             , fitting_done(n_atoms))
   basis_done   = .false.
   fitting_done = .false.
   min_exp      = 100000.0D0
   n_orig       = 0
   n_aux        = 0

   open(unit = file_uid, file= basis_file, iostat = file_iostat)
   read(file_uid,'(A8)') start_str
   do while (start_str .ne. 'endbasis')
      ! Starts reading the basis set for an atom
      ! Reads contraction scheme. The value for p/d/f should not be repeated
      ! 3/6/10 times, just set them once. Also reads angular momentum for each
      ! of the contractions.
      read(file_uid,*) atom, nraw, ncon
      read(file_uid,*) n_cont_func(1:ncon)
      read(file_uid,*) ang_mom(1:ncon)
      do icount = 1, nraw
         read(file_uid,*) expo_temp(icount), coef_temp(icount)
         if ((expo_temp(icount) .lt. min_exp) .and. any(atom_Z == atom)) &
                                                    min_exp = expo_temp(icount)
      enddo

      do iatom = 1, n_atoms
         if (atom_Z(iatom) .eq. atom .and. (.not. basis_done(iatom))) then
            basis_done(iatom)  = .true.
            min_atm_exp(iatom) = min_exp

            ! These are used for atoms that are near to each other.
            nns(iatom) = 0
            nnp(iatom) = 0
            nnd(iatom) = 0

            do icont = 1, ncon
               select case (ang_mom(icont))
               case (0)
                  nns(iatom) = nns(iatom) + ANG_DEG(ang_mom(icont))
               case (1)
                  nnp(iatom) = nnp(iatom) + ANG_DEG(ang_mom(icont))
               case (2)
                  nnd(iatom) = nnd(iatom) + ANG_DEG(ang_mom(icont))
               case default
               end select
            enddo

            ! Stores the exponents and coefficients, and normalizes them if
            ! necessary.
            index  = 0
            do icont = 1, ncon
               nshell(ang_mom(icont)) = nshell(ang_mom(icont)) + &
                                        ANG_DEG(ang_mom(icont))
               do l2 = 1, ANG_DEG(ang_mom(icont))
                  n_orig = n_orig +1

                  if (normalize) then
                     do icount =1, n_cont_func(icont)
                        index = index +1
                        select case (ang_mom(icont))
                        case (0)
                           craw(n_orig, icount) = coef(n_orig, icount)
                           coef(n_orig, icount) = dsqrt( dsqrt(8.0D0 * &
                                                  (expo_temp(index)) ** 3 ) / &
                                                  PI32) * coef_temp(index)
                           expo(n_orig, icount) = expo_temp(index)
                        case (1)
                           craw(n_orig, icount) = coef(n_orig, icount)
                           coef(n_orig, icount) = dsqrt( dsqrt(8.0D0 * &
                                                  (expo_temp(index)) ** 3 ) * &
                                                  4.0D0 * expo_temp(index) /  &
                                                  PI32) * coef_temp(index)
                           expo(n_orig, icount) = expo_temp(index)
                        case (2)
                           craw(n_orig, icount) = coef(n_orig, icount)
                           coef(n_orig, icount) = dsqrt( dsqrt(8.0D0 * &
                                                  (expo_temp(index)) ** 3 ) *  &
                                                  16.0D0 * expo_temp(index)**2/&
                                                  PI32) * coef_temp(index)
                           expo(n_orig, icount) = expo_temp(index)
                        case default
                           iostatus = 1
                           write(*,'(A)') "  ERROR: Basis set contains "    , &
                                          "usupported angular momenta. LIO ", &
                                          "uses only s, p and d-type orbitals."
                           return
                        end select
                     enddo
                  else
                     do icount =1, n_cont_func(icont)
                        index = index +1
                        coef(n_orig, icount) = coef_temp(index)
                        expo(n_orig, icount) = expo_temp(index)
                     enddo
                  endif

                  ! Repeats the index for p, d and f
                  if (l2 .ne. ANG_DEG(ang_mom(icont))) index = index - &
                                                              n_cont_func(icont)

                  atm_of_func(n_orig) = iatom
                  n_cont(n_orig)       = n_cont_func(icont)
                  ang_mom_f(n_orig)   = ang_mom(icont)
               enddo
            enddo
         endif
      enddo

      ! Starts reading the auxiliary basis for an atom. Same rules as before
      ! are applied here.
      read(file_uid,*) atom, nraw, ncon
      read(file_uid,*) n_cont_func(1:ncon)
      read(file_uid,*) ang_mom(1:ncon)
      do icount = 1, nraw
         read(file_uid,*) expo_temp(icount), coef_temp(icount)
      enddo

      do iatom = 1, n_atoms
         if (atom_Z(iatom) .eq. atom) then
         fitting_done(iatom) = .true.

         index = 0
         do icont = 1, ncon
            nshelld(ang_mom(icont)) = nshelld(ang_mom(icont)) + &
                                      ANG_DEG(ang_mom(icont))
            do l2 = 1, ANG_DEG(ang_mom(icont))
               n_aux = n_aux +1

               if (normalize) then
                  do icount =1, n_cont_func(icont)
                     index = index +1
                     select case (ang_mom(icont))
                     case (0)
                        coefd(n_aux, icount) = dsqrt( dsqrt(8.0D0 * &
                                               (expo_temp(index)) ** 3 ) / &
                                               PI32) * coef_temp(index)
                        expod(n_aux, icount) = expo_temp(index)
                     case (1)
                        coefd(n_aux, icount) = dsqrt( dsqrt(8.0D0 * &
                                               (expo_temp(index)) ** 3 ) * &
                                               4.0D0 * expo_temp(index) /  &
                                               PI32) * coef_temp(index)
                        expod(n_aux, icount) = expo_temp(index)
                     case (2)
                        coefd(n_aux, icount) = dsqrt( dsqrt(8.0D0 * &
                                               (expo_temp(index)) ** 3 ) *   &
                                               16.0D0 * expo_temp(index) **2/&
                                               PI32) * coef_temp(index)
                        expod(n_aux, icount) = expo_temp(index)
                     case default
                        iostatus = 1
                        write(*,'(A)') "  ERROR: Basis set contains "    , &
                                       "usupported angular momenta. LIO ", &
                                       "uses only s, p and d-type orbitals."
                        return
                     end select
                  enddo
               else
                  do icount =1, n_cont_func(icont)
                     index = index +1
                     coefd(n_aux, icount) = coef_temp(index)
                     expod(n_aux, icount) = expo_temp(index)
                  enddo
               endif

               if (l2 .ne. ANG_DEG(ang_mom(icont))) then
                  index = index - n_cont_func(icont)
               endif

               atm_of_funcd(n_aux) = iatom
               n_contd(n_aux)      = n_cont_func(icont)
               ang_mom_fd(n_aux)   = ang_mom(icont)
            enddo
         enddo
       endif
   enddo
   read(file_uid,'(A8)') start_str
   enddo
end subroutine read_basis_external

subroutine read_basis_internal(basis_file, fitting_file, n_funcs, n_fits,     &
                               n_atoms, normalize, max_fun_per_atom,          &
                               max_con_per_atom, atom_Z, craw, coef, expo,    &
                               coefd, expod, n_cont, n_contd, ang_mom_f,      &
                               ang_mom_fd, atm_of_func, atm_of_funcd,         &
                               min_atm_exp, nns, nnp, nnd, nShell, nShelld,   &
                               iostatus)
   use basis_data   , only: ANG_DEG
   use constants_mod, only: PI32
   implicit none
   integer         , intent(in)  :: max_con_per_atom, max_fun_per_atom, &
                                    n_atoms, atom_Z(n_atoms), n_funcs, n_fits
   logical         , intent(in)  :: normalize
   character(len=*), intent(in)  :: basis_file, fitting_file
   integer         , intent(out) :: atm_of_func(n_funcs), atm_of_funcd(n_fits),&
                                    nns(n_atoms), nnp(n_atoms), nnd(n_atoms),  &
                                    iostatus, n_cont(n_funcs), n_contd(n_fits),&
                                    ang_mom_f(n_funcs),  ang_mom_fd(n_fits),   &
                                    nShell(0:3), nShelld(0:3)
   double precision, intent(out) :: coef(n_funcs,max_con_per_atom), &
                                    craw(n_funcs,max_con_per_atom), &
                                    expo(n_funcs,max_con_per_atom), &
                                    coefd(n_fits,max_con_per_atom), &
                                    expod(n_fits,max_con_per_atom), &
                                    min_atm_exp(n_atoms)

   integer            :: file_iostat, file_uid = TMP_OPEN_UID
   integer            :: iatom, nraw, ncon, atom, icont, icount, l2, index,    &
                         n_orig, n_aux
   double precision   :: min_exp
   character(len=20)  :: start_str
   character(len=100) :: line_read

   logical         , allocatable :: basis_done(:), fitting_done(:)
   integer         , allocatable :: n_cont_func(:), ang_mom(:)
   double precision, allocatable :: expo_temp(:), coef_temp(:)

   allocate(n_cont_func(max_con_per_atom +1), ang_mom(max_con_per_atom +1), &
            expo_temp(max_fun_per_atom)     , coef_temp(max_fun_per_atom) , &
            basis_done(n_atoms)             , fitting_done(n_atoms))
   basis_done   = .false.
   fitting_done = .false.
   min_exp      = 100000.0D0
   n_orig       = 0
   n_aux        = 0

   open(unit = file_uid, file= basis_file, iostat = file_iostat)
   ! Skips # characters in file.
   line_read = ""
   do while (line_read == "")
      read(file_uid,*) line_read
      read(line_read, '(A1)') start_str
      if (trim(start_str) == "#") line_read = ""
   enddo
   read(line_read, '(A8)') start_str

   do while (start_str .ne. 'endbasis')
      ! Starts reading the basis set for an atom
      ! Reads contraction scheme. The value for p/d/f should not be repeated
      ! 3/6/10 times, just set them once. Also reads angular momentum for each
      ! of the contractions.
      read(file_uid,*) atom, nraw, ncon
      read(file_uid,*) n_cont_func(1:ncon)
      read(file_uid,*) ang_mom(1:ncon)
      do icount = 1, nraw
         read(file_uid,*) expo_temp(icount), coef_temp(icount)
         if ((expo_temp(icount) .lt. min_exp) .and. any(atom_Z == atom)) &
                                                    min_exp = expo_temp(icount)
      enddo

      do iatom = 1, n_atoms
         if (atom_Z(iatom) .eq. atom .and. (.not. basis_done(iatom))) then
            basis_done(iatom)  = .true.
            min_atm_exp(iatom) = min_exp

            ! These are used for atoms that are near to each other.
            nns(iatom) = 0
            nnp(iatom) = 0
            nnd(iatom) = 0

            do icont = 1, ncon
               select case (ang_mom(icont))
               case (0)
                  nns(iatom) = nns(iatom) + ANG_DEG(ang_mom(icont))
               case (1)
                  nnp(iatom) = nnp(iatom) + ANG_DEG(ang_mom(icont))
               case (2)
                  nnd(iatom) = nnd(iatom) + ANG_DEG(ang_mom(icont))
               case default
               end select
            enddo

            ! Stores the exponents and coefficients, and normalizes them if
            ! necessary.
            index  = 0
            do icont = 1, ncon
               nshell(ang_mom(icont)) = nshell(ang_mom(icont)) + &
                                        ANG_DEG(ang_mom(icont))
               do l2 = 1, ANG_DEG(ang_mom(icont))
                  n_orig = n_orig +1

                  if (normalize) then
                     do icount =1, n_cont_func(icont)
                        index = index +1
                        select case (ang_mom(icont))
                        case (0)
                           craw(n_orig, icount) = coef(n_orig, icount)
                           coef(n_orig, icount) = dsqrt( dsqrt(8.0D0 * &
                                                  (expo_temp(index)) ** 3 ) / &
                                                  PI32) * coef_temp(index)
                           expo(n_orig, icount) = expo_temp(index)
                        case (1)
                           craw(n_orig, icount) = coef(n_orig, icount)
                           coef(n_orig, icount) = dsqrt( dsqrt(8.0D0 * &
                                                  (expo_temp(index)) ** 3 ) * &
                                                  4.0D0 * expo_temp(index) /  &
                                                  PI32) * coef_temp(index)
                           expo(n_orig, icount) = expo_temp(index)
                        case (2)
                           craw(n_orig, icount) = coef(n_orig, icount)
                           coef(n_orig, icount) = dsqrt( dsqrt(8.0D0 * &
                                                  (expo_temp(index)) ** 3 ) *  &
                                                  16.0D0 * expo_temp(index)**2/&
                                                  PI32) * coef_temp(index)
                           expo(n_orig, icount) = expo_temp(index)
                        case default
                           iostatus = 1
                           write(*,'(A)') "  ERROR: Basis set contains "    , &
                                          "usupported angular momenta. LIO ", &
                                          "uses only s, p and d-type orbitals."
                           return
                        end select
                     enddo
                  else
                     do icount =1, n_cont_func(icont)
                        index = index +1
                        coef(n_orig, icount) = coef_temp(index)
                        expo(n_orig, icount) = expo_temp(index)
                     enddo
                  endif

                  ! Repeats the index for p, d and f
                  if (l2 .ne. ANG_DEG(ang_mom(icont))) index = index - &
                                                              n_cont_func(icont)

                  atm_of_func(n_orig) = iatom
                  n_cont(n_orig)      = n_cont_func(icont)
                  ang_mom_f(n_orig)   = ang_mom(icont)
               enddo
            enddo
         endif
      enddo

      ! Skips # characters in file.
      line_read = ""
      do while (line_read == "")
         read(file_uid,*) line_read
         read(line_read, '(A1)') start_str
         if (trim(start_str) == "#") line_read = ""
      enddo
      read(line_read, '(A8)') start_str
   enddo
   close(file_uid)

   open(unit = file_uid, file= fitting_file, iostat = file_iostat)
   ! Skips # characters in file.
   line_read = ""
   do while (line_read == "")
      read(file_uid,*) line_read
      read(line_read, '(A1)') start_str
      if (trim(start_str) == "#") line_read = ""
   enddo
   read(line_read, '(A8)') start_str

   do while (start_str .ne. 'endbasis')
      ! Starts reading the auxiliary basis for an atom. Same rules as before
      ! are applied here.
      read(file_uid,*) atom, nraw, ncon
      read(file_uid,*) n_cont_func(1:ncon)
      read(file_uid,*) ang_mom(1:ncon)
      do icount = 1, nraw
         read(file_uid,*) expo_temp(icount), coef_temp(icount)
      enddo

      do iatom = 1, n_atoms
         if (atom_Z(iatom) .eq. atom) then
            fitting_done(iatom) = .true.

            index = 0
            do icont = 1, ncon
               nshelld(ang_mom(icont)) = nshelld(ang_mom(icont)) + &
                                         ANG_DEG(ang_mom(icont))
               do l2 = 1, ANG_DEG(ang_mom(icont))
                  n_aux = n_aux +1

                  if (normalize) then
                     do icount =1, n_cont_func(icont)
                        index = index +1
                        select case (ang_mom(icont))
                        case (0)
                           coefd(n_aux, icount) = dsqrt( dsqrt(8.0D0 * &
                                                  (expo_temp(index)) ** 3 ) / &
                                                  PI32) * coef_temp(index)
                           expod(n_aux, icount) = expo_temp(index)
                        case (1)
                           coefd(n_aux, icount) = dsqrt( dsqrt(8.0D0 * &
                                                  (expo_temp(index)) ** 3 ) * &
                                                  4.0D0 * expo_temp(index) /  &
                                                  PI32) * coef_temp(index)
                           expod(n_aux, icount) = expo_temp(index)
                        case (2)
                           coefd(n_aux, icount) = dsqrt( dsqrt(8.0D0 * &
                                                  (expo_temp(index)) ** 3 ) *  &
                                                  16.0D0 * expo_temp(index)**2/&
                                                  PI32) * coef_temp(index)
                           expod(n_aux, icount) = expo_temp(index)
                        case default
                           iostatus = 1
                           write(*,'(A)') "  ERROR: Basis set contains "    , &
                                          "usupported angular momenta. LIO ", &
                                          "uses only s, p and d-type orbitals."
                           return
                        end select
                     enddo
                  else
                     do icount =1, n_cont_func(icont)
                        index = index +1
                        coefd(n_aux, icount) = coef_temp(index)
                        expod(n_aux, icount) = expo_temp(index)
                     enddo
                  endif

                  if (l2 .ne. ANG_DEG(ang_mom(icont))) then
                     index = index - n_cont_func(icont)
                  endif

                  atm_of_funcd(n_aux) = iatom
                  n_contd(n_aux)      = n_cont_func(icont)
                  ang_mom_fd(n_aux)   = ang_mom(icont)
               enddo
            enddo
          endif
      enddo

      ! Skips # characters in file.
      line_read = ""
      do while (line_read == "")
         read(file_uid,*) line_read
         read(line_read, '(A1)') start_str
         if (trim(start_str) == "#") line_read = ""
      enddo
      read(line_read, '(A8)') start_str
   enddo

end subroutine read_basis_internal

subroutine reorder_basis(expon, coeff, atom_of_funct, n_cont, mixed_index, &
                         basis_size, max_cont, l_of_funct, n_shell, craw)
   implicit none
   integer         , intent(in)    :: basis_size, max_cont, n_shell(0:3), &
                                      l_of_funct(basis_size)
   integer         , intent(inout) :: atom_of_funct(basis_size),   &
                                      n_cont(basis_size)      ,    &
                                      mixed_index(basis_size)
   double precision, intent(inout) :: expon(basis_size, max_cont), &
                                      coeff(basis_size, max_cont)
   double precision, intent(inout), optional :: craw(basis_size, max_cont)

   double precision, allocatable :: expo_t(:,:), coef_t(:,:), craw_t(:,:)
   integer         , allocatable :: atom_of_funct_t(:), n_cont_t(:)
   integer :: ifunct, s_index, p_index, d_index

   allocate(expo_t(basis_size, max_cont), coef_t(basis_size, max_cont), &
            atom_of_funct_t(basis_size) , n_cont_t(basis_size)        , &
            craw_t(basis_size, max_cont) )

   s_index = 1
   p_index = 1 + n_shell(0)
   d_index = 1 + n_shell(0) + n_shell(1)

   do ifunct = 1, basis_size
      select case (l_of_funct(ifunct))
      case (0) ! s function
         atom_of_funct_t(s_index) = atom_of_funct(ifunct)
         mixed_index(s_index)     = ifunct
         n_cont_t(s_index)        = n_cont(ifunct)
         expo_t(s_index,:)        = expon(ifunct,:)
         coef_t(s_index,:)        = coeff(ifunct,:)
         if (present(craw)) craw_t(s_index,:) = craw(ifunct,:)

         s_index = s_index +1
      case (1) ! p functions
         atom_of_funct_t(p_index) = atom_of_funct(ifunct)
         mixed_index(p_index)     = ifunct
         n_cont_t(p_index)        = n_cont(ifunct)
         expo_t(p_index,:)        = expon(ifunct,:)
         coef_t(p_index,:)        = coeff(ifunct,:)
         if (present(craw)) craw_t(p_index,:) = craw(ifunct,:)

         p_index = p_index +1
      case (2) ! d functions
         atom_of_funct_t(d_index) = atom_of_funct(ifunct)
         mixed_index(d_index)     = ifunct
         n_cont_t(d_index)        = n_cont(ifunct)
         expo_t(d_index,:)        = expon(ifunct,:)
         coef_t(d_index,:)        = coeff(ifunct,:)
         if (present(craw)) craw_t(d_index,:) = craw(ifunct,:)

         d_index = d_index +1
      case default
      end select
   enddo

   n_cont        = n_cont_t
   atom_of_funct = atom_of_funct_t
   expon         = expo_t
   coeff         = coef_t
   if (present(craw)) craw = craw_t
   deallocate(expo_t, coef_t, atom_of_funct_t, n_cont_t, craw_t)
end subroutine reorder_basis
end module basis_subs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
