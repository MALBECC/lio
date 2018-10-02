!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% BASIS.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains two modules: basis_data, containing basis set function    !
! data, and basis_subs, containing the following subroutines:                  !
!                                                                              !
! · basis_set_size(): Prereads data and sets array sizes                       !
! · basis_read()    : Reads input basis functions.                             !
! · basis_init()    : Initializes basis data.                                  !
! · basis_deinit()  : Deinitializes basis data.                                !
! · funct()         : Boys function.                                           !
!                                                                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

! TODO: INVERT C and A matrices, since loops are faster that way. In Ehrenfest
!       these are already inverted.

module basis_data
   implicit none
   ! Namelist inputfiles
   ! int_basis : Use LIO internal basis (located in /dat)
   ! rMax      : Maximum exponent for double-precision integrals.
   ! rMaxs     : Maximum exponent for single-precision integrals.
   ! norm      : Normalize integrals (deprecated).
   logical          :: int_basis = .true.
   double precision :: rMax      = 16.0D0
   double precision :: rMaxs     =  5.0D0
   logical          :: norm      = .true.

   ! Single variables
   ! M     : number of basis functions.
   ! Md    : number of fitting functions.
   ! kknums: number of single-precision two-center integrals.
   ! kknumd: number of double-precision two-center integrals.
   integer          :: M      = 0
   integer          :: Md     = 0
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

   ! Degeneracy for each angular momentum
   integer         , parameter :: ANG_DEG(0:3) = (/1, 3, 6, 10/)
contains
end module basis_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module basis_subs
   implicit none
   integer, parameter :: TMP_OPEN_UID = 9999
contains

subroutine basis_init(basis_name, fitting_name, n_atoms, atom_Z, out_stat)
   !use basis_data, only: M, Md, int_basis, Nuc, Nucd, nCont, nContd, a, c, ad, &
   use garcha_mod, only: M, Md, int_basis, Nuc, Nucd, nCont, nContd, a, c, ad, &
                         cd, atmin, nns, nnp, nnd, nshell, nshelld, norm
   use basis_data, only: ang_mom, ang_momd, max_f_per_atom, max_c_per_atom
   implicit none
   integer         , intent(in)    :: n_atoms, atom_Z(n_atoms)
   integer         , intent(out)   :: out_stat
   character(len=*), intent(inout) :: basis_name, fitting_name

   integer :: icount, iostat
   integer, allocatable :: atom_count(:)
   logical, allocatable :: atom_basis_chk(:), atom_fitting_chk(:)

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
   !print*, "DBG ", "M = ", M, "Md = ", Md, "max_c = ", max_c_per_atom, &
   !        "max_f = ", max_f_per_atom
   call check_basis(atom_basis_chk, atom_fitting_chk, n_atoms, atom_Z, iostat)
   if (iostat .gt. 0) then
      out_stat = 2
      return
   endif

   allocate(c(M, max_c_per_atom), a(M, max_c_per_atom), cd(Md, max_c_per_atom),&
            ad(Md, max_c_per_atom), nCont(M), nContd(Md), ang_mom(M),        &
            ang_momd(Md), Nuc(M), Nucd(M), atmin(n_atoms), nns(n_atoms),     &
            nnp(n_atoms), nnd(n_atoms))
   if (int_basis) then
      call read_basis_internal(basis_name, fitting_name, M, Md, n_atoms, norm, &
                               max_f_per_atom, max_c_per_atom, atom_Z, c, a,   &
                               cd, ad, nCont, nContd, ang_mom, ang_momd, Nuc,  &
                               Nucd, atmin, nns, nnp, nnd, nshell, nshelld,    &
                               iostat)
   else
      call read_basis_external(basis_name, M,Md, n_atoms, norm, max_f_per_atom,&
                               max_c_per_atom, atom_Z, c, a, cd, ad, nCont,    &
                               nContd, ang_mom, ang_momd, Nuc, Nucd, atmin,    &
                               nns, nnp, nnd, nshell, nshelld, iostat)
   endif
   if (iostat .gt. 0) then
      out_stat = 3
      return
   endif

   deallocate(atom_count, atom_basis_chk, atom_fitting_chk, ang_mom, ang_momd)
end subroutine basis_init

subroutine basis_setup_ehren()
   use basis_data, only: a_ehren, c_ehren, ang_mom_ehren, max_c_per_atom, &
                         max_f_per_atom
   use garcha_mod, only: a, c, nShell, M

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
   enddo
end subroutine basis_setup_ehren

subroutine basis_set_size(basis_size, aux_size, max_f_per_atom, max_c_per_atom,&
                          atom_bas_done, atom_fit_done, basis_file,            &
                          fitting_file, atom_count, atom_Z, n_atoms,           &
                          use_internal, iostatus)
   use basis_data, only: ANG_DEG
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
   integer              :: file_iostat, icount
   integer              :: iatom, nraw, ncon, l_of_func
   character(len=20)    :: start_str
   character(len=100)   :: lio_dir, line_read

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

         ! Only needs to read angular momenta for each contraction, the rest is
         ! skipped.
         read(file_uid,*)
         do icount = 1, ncon
            read(file_uid, *) l_of_func
            if (any(atom_Z == iatom)) then
               basis_size = basis_size + ANG_DEG(l_of_func) * atom_count(iatom)
            endif
         enddo
         do icount = 1, nraw
            read(file_uid,*)
         enddo
         atom_bas_done(iatom) = .true.

         ! Reads auxiliary basis set for an atom doing the same as before.
         read(file_uid,*) iatom, nraw, ncon
         if (max_f_per_atom .lt. nraw) max_f_per_atom = nraw
         if (max_c_per_atom .lt. ncon) max_c_per_atom = ncon

         read(file_uid,*)
         do icount = 1, ncon
            read(file_uid,*) l_of_func
            if (any(atom_Z == iatom)) then
               aux_size = aux_size + ANG_DEG(l_of_func) * atom_count(iatom)
            endif
         enddo
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

         ! Only needs to read angular momenta for each contraction, the rest is
         ! skipped.
         read(file_uid,*)
         do icount = 1, ncon
            read(file_uid,*) l_of_func
            if (any(atom_Z == iatom)) then
               basis_size = basis_size + ANG_DEG(l_of_func) * atom_count(iatom)
            endif
         enddo
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
      open(unit = file_uid, file= basis_file, iostat = file_iostat)
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

         ! Only needs to read angular momenta for each contraction, the rest is
         ! skipped.
         read(file_uid,*)
         do icount = 1, ncon
            read(file_uid,*) l_of_func
            if (any(atom_Z == iatom)) then
               aux_size = aux_size + ANG_DEG(l_of_func) * atom_count(iatom)
            endif
         enddo
         read(file_uid,*)
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
end subroutine basis_set_size

subroutine check_basis(atom_bas_done, atom_fit_done, n_atoms, atom_Z, iostatus)
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
                               coef, expo, coefd, expod, n_cont, n_contd,      &
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
   integer         , intent(out) :: atm_of_func(n_funcs),atm_of_funcd(n_funcs),&
                                    nns(n_atoms), nnp(n_atoms), nnd(n_atoms),  &
                                    iostatus, n_cont(n_funcs), n_contd(n_fits),&
                                    ang_mom_f(n_funcs),  ang_mom_fd(n_fits),   &
                                    nShell(0:3), nShelld(0:3)
   double precision, intent(out) :: coef(n_funcs,max_con_per_atom), &
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
   read(file_uid,'(A8)') start_str
   do while (start_str .ne. 'endbasis')
      ! Starts reading the basis set for an atom
      ! Reads contraction scheme. The value for p/d/f should not be repeated
      ! 3/6/10 times, just set them once. Also reads angular momentum for each
      ! of the contractions.
      read(file_uid,*) atom, nraw, ncon
      read(file_uid,*) n_cont_func
      read(file_uid,*) ang_mom
      do icount = 1, nraw
         read(file_uid,*) expo_temp(icount), coef_temp(icount)
         if (expo_temp(icount) .lt. min_exp) min_exp = expo_temp(icount)
      enddo

      do iatom = 1, n_atoms
         if (atom_Z(iatom) .eq. atom .and. (.not. basis_done(iatom))) then
            basis_done(iatom)   = .true.
            min_atm_exp(iatom) = min_exp

            ! These are used for atoms that are near to each other.
            nns(iatom) = 0
            nnp(iatom) = 0
            nnd(iatom) = 0

            do icont = 1, ncon
               select case (ang_mom(icont))
               case (0)
                  nns(iatom) = nns(iatom) + ANG_DEG(ang_mom(icont))
               case (2)
                  nnp(iatom) = nnp(iatom) + ANG_DEG(ang_mom(icont))
               case (3)
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
                           coef(n_orig, icount) = dsqrt( dsqrt(8.0D0 * &
                                                  (expo_temp(index)) ** 3 ) / &
                                                  PI32) * coef_temp(index)
                           expo(n_orig, icount) = expo_temp(index)
                        case (1)
                           coef(n_orig, icount) = dsqrt( dsqrt(8.0D0 * &
                                                  (expo_temp(index)) ** 3 ) * &
                                                  4.0D0 * expo_temp(index) /  &
                                                  PI32) * coef_temp(index)
                           expo(n_orig, icount) = expo_temp(index)
                        case (2)
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
      read(file_uid,*) n_cont_func
      read(file_uid,*) ang_mom
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
                               max_con_per_atom, atom_Z, coef, expo, coefd,   &
                               expod, n_cont, n_contd, ang_mom_f, ang_mom_fd, &
                               atm_of_func, atm_of_funcd, min_atm_exp, nns,   &
                               nnp, nnd, nShell, nShelld, iostatus)
   use basis_data   , only: ANG_DEG
   use constants_mod, only: PI32
   implicit none
   integer         , intent(in)  :: max_con_per_atom, max_fun_per_atom, &
                                    n_atoms, atom_Z(n_atoms), n_funcs, n_fits
   logical         , intent(in)  :: normalize
   character(len=*), intent(in)  :: basis_file, fitting_file
   integer         , intent(out) :: atm_of_func(n_funcs),atm_of_funcd(n_funcs),&
                                    nns(n_atoms), nnp(n_atoms), nnd(n_atoms),  &
                                    iostatus, n_cont(n_funcs), n_contd(n_fits),&
                                    ang_mom_f(n_funcs),  ang_mom_fd(n_fits),   &
                                    nShell(0:3), nShelld(0:3)
   double precision, intent(out) :: coef(n_funcs,max_con_per_atom), &
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
      read(file_uid,*) n_cont_func
      read(file_uid,*) ang_mom
      do icount = 1, nraw
         read(file_uid,*) expo_temp(icount), coef_temp(icount)
         if (expo_temp(icount) .lt. min_exp) min_exp = expo_temp(icount)
      enddo

      do iatom = 1, n_atoms
         if (atom_Z(iatom) .eq. atom .and. (.not. basis_done(iatom))) then
            basis_done(iatom)   = .true.
            min_atm_exp(iatom) = min_exp

            ! These are used for atoms that are near to each other.
            nns(iatom) = 0
            nnp(iatom) = 0
            nnd(iatom) = 0

            do icont = 1, ncon
               select case (ang_mom(icont))
               case (0)
                  nns(iatom) = nns(iatom) + ANG_DEG(ang_mom(icont))
               case (2)
                  nnp(iatom) = nnp(iatom) + ANG_DEG(ang_mom(icont))
               case (3)
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
                           coef(n_orig, icount) = dsqrt( dsqrt(8.0D0 * &
                                                  (expo_temp(index)) ** 3 ) / &
                                                  PI32) * coef_temp(index)
                           expo(n_orig, icount) = expo_temp(index)
                        case (1)
                           coef(n_orig, icount) = dsqrt( dsqrt(8.0D0 * &
                                                  (expo_temp(index)) ** 3 ) * &
                                                  4.0D0 * expo_temp(index) /  &
                                                  PI32) * coef_temp(index)
                           expo(n_orig, icount) = expo_temp(index)
                        case (2)
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
      read(file_uid,*) n_cont_func
      read(file_uid,*) ang_mom
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




!###############################################################################
! TAKE THIS TO OTHER MODULE
subroutine atom_name(atom_Z, symb)
 ! Takes atomic number Z and translates it to its name.
 implicit none
 integer         , intent(in)  :: atom_Z
 character(LEN=3), intent(out) :: symb

 character(LEN=3) :: name(118)
 name = (/'H  ','HE ','LI ','BE ','B  ','C  ','N  ','O  ','F  ','NE ','NA ', &
          'MG ','AL ','SI ','P  ','S  ','CL ','AR ','K  ','CA ','SC ','TI ', &
          'V  ','CR ','MN ','FE ','CO ','NI ','CU ','ZN ','GA ','GE ','AS ', &
          'SE ','BR ','KR ','RB ','SR ','Y  ','ZR ','NB ','MO ','TC ','RU ', &
          'RH ','PD ','AG ','CD ','IN ','SN ','SB ','TE ','I  ','XE ','CS ', &
          'BA ','LA ','CE ','PR ','ND ','PM ','SM ','EU ','GD ','TB ','DY ', &
          'HO ','ER ','TM ','YB ','LU ','HF ','TA ','W  ','RE ','OS ','IR ', &
          'PT ','AU ','HG ','TL ','PB ','BI ','PO ','AT ','RN ','FR ','RA ', &
          'AC ','TH ','PA ','U  ','NP ','PU ','AM ','CM ','BK ','CF ','ES ', &
          'FM ','MD ','NO ','LR ','RF ','DB ','SG ','BH ','HS ','MT ','DS ', &
          'UUU','UUB','UUT','UUQ','UUP','UUH','UUS','UUO'/)
 symb = name(atom_Z)

end subroutine atom_name

end module basis_subs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
