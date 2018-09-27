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
module basis_data
   implicit none
   ! Namelist inputfiles
   logical          :: int_basis = .true.
   double precision :: rmax      = 16.0D0
   double precision :: rmaxs     =  5.0D0
   logical          :: NORM      = .true.

   ! Single variables
   integer          :: M      = 0
   integer          :: Md     = 0
   integer          :: kknums = 0
   integer          :: kknumd = 0

   ! Preset arrays
   integer :: nShell(0:4)
   integer :: nShelld(0:4)

   ! Allocatable arrays
   integer         , allocatable :: Nuc(:)
   integer         , allocatable :: Nucd(:)
   integer         , allocatable :: ncont(:)
   integer         , allocatable :: ncontd(:)
   integer         , allocatable :: kkind(:)
   integer         , allocatable :: kkinds(:)
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
   double precision, allocatable :: cool(:)
   real            , allocatable :: cools(:)

   ! Degeneracy for each angular momentum
   double precision :: ang_deg(0:4) = (/1, 3, 6, 10/)
contains
end module basis_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module basis_subs
contains

subroutine basis_init(basis_name, fitting_name, n_atoms, atom_Z, out_stat)
   use basis_data, only: M, Md, int_basis
   implicit none
   integer         , intent(in)  :: n_atoms, atom_Z(n_atoms)
   character(len=*), intent(in)  :: basis_name, fitting_name
   integer         , intent(out) :: out_stat

   integer :: iostat, max_f_per_atom
   integer, allocatable :: atom_count(:)
   logical, allocatable :: atom_basis_chk(:), atom_fitting_chk(:)

   allocate(atom_count(0:120), atom_basis_chk(0:120), atom_fitting_chk(0:120))
   ! Precalculates the amount of atoms of a certain type, in order to
   ! correctly asses the sizes of M and Md
   atom_count = 0
   do icount = 1, n_atoms
      atom_count(Iz(icount)) = atom_count(Iz(icount)) + 1
   enddo

   call basis_set_size(M, Md, max_f_per_atom, max_c_per_atom, atom_basis_chk,  &
                       atom_fitting_chk, basis_name, fitting_name, atom_count, &
                       atom_Z, n_atoms, use_internal, iostat)
   if (iostat .gt. 0) then
      out_stat = 1
      return
   endif
   !print*, "DBG ", "M = ", M, "Md = ", Md, "max_c = ", max_c_per_atom, &
   !        "max_f = ", max_f_per_atom

   call check_basis(atom_basis_chk, atom_fitting_chk, iostat, n_atoms, atom_Z)
   if (iostat .gt. 0) then
      out_stat = 2
      return
   endif

   if (int_basis) then
      call read_basis_internal()
   else
      call read_basis_external()
   endif

   ! allocatear ncf(max_c_per_atom+1), lt(max_c_per_atom+1)
   ! allocatear at(max_func), ct(max_func)

   deallocate(atom_count, atom_basis_chk, atom_fitting_chk)
end subroutine basis_init

subroutine basis_set_size(basis_size, aux_size, max_f_per_atom, max_c_per_atom,&
                          basis_file, fitting_file, atom_count, atom_Z,        &
                          n_atoms, use_internal, iostatus)
   use basis_data, only: ang_deg
   implicit none
   integer         , intent(in)  :: n_atoms, max_f_per_atom, atom_Z(n_atoms), &
                                    atom_count(0:120)
   logical         , intent(in)  :: use_internal
   character(len=*), intent(in)  :: basis_file, fitting_file
   integer         , intent(out) :: basis_size, aux_size, max_f_per_atom, &
                                    max_c_per_atom, iostatus
   logical         , intent(out) :: atom_bas_done(0:120), atom_fit_done(0:120)

   logical              :: file_exists
   integer              :: file_iostat, icount
   integer              :: iatom, nraw, ncon, ang_mom
   character(len=20)    :: start_str
   character(len=*)     :: lio_dir, line_read

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
      open(unit = 9999, file= basis_file, iostat = file_iostat)
      read(9999, '(A8)') start_str

      do while (start_str .ne. "endbasis")
         read(9999,*) iatom, nraw, ncon
         if (max_f_per_atom .lt. nraw) max_f_per_atom = nraw
         if (max_c_per_atom .lt. ncon) max_c_per_atom = ncon

         ! Only needs to read angular momenta for each contraction, the rest is
         ! skipped.
         do icount = 1, ncon
            read(9999,*)
         enddo
         do icount = 1, ncon
            read(9999,*) ang_mom
            if (any(atom_Z == iatom)) then
               basis_size = basis_size + ang_deg(ang_mom) * atom_count(iatom)
            endif
         enddo
         do icount = 1, nraw
            read(9999,*)
         enddo
         atom_bas_done(iatom) = .true.

         ! Reads auxiliary basis set for an atom doing the same as before.
         read(9999,*) iatom, nraw, ncon
         if (max_f_per_atom .lt. nraw) max_f_per_atom = nraw
         if (max_c_per_atom .lt. ncon) max_c_per_atom = ncon

         do icount = 1, ncon
            read(9999,*)
         enddo
         do icount = 1, ncon
            read(9999,*) ang_mom
            if (any(atom_Z == iatom)) then
               aux_size = aux_size + ang_deg(ang_mom) * atom_count(iatom)
            endif
         enddo
         do icount = 1, nraw
            read(9999,*)
         enddo
         atom_fit_done(iatom) = .true.

         read(9999, '(A8)') start_str
      enddo
      close(9999)
   else
      ! Reads from LIO internal basis set files.
      call getenv("LIOHOME", lio_dir)
      if (lio_dir = "") then
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
      open(unit = 9999, file= basis_file, iostat = file_iostat)
      line_read = ""
      do while (line_read == "")
         read(9999,*) line_read
         read(line_read, '(A1)') start_str
         if (start_str = "#") line_read = ""
      enddo
      read(line_read, '(A8)') start_str

      do while (start_str .ne. "endbasis")
         read(9999,*) iatom, nraw, ncon
         if (max_f_per_atom .lt. nraw) max_f_per_atom = nraw
         if (max_c_per_atom .lt. ncon) max_c_per_atom = ncon

         ! Only needs to read angular momenta for each contraction, the rest is
         ! skipped.
         do icount = 1, ncon
            read(9999,*)
         enddo
         do icount = 1, ncon
            read(9999,*) ang_mom
            if (any(atom_Z == iatom)) then
               basis_size = basis_size + ang_deg(ang_mom) * atom_count(iatom)
            endif
         enddo
         do icount = 1, nraw
            read(9999,*)
         enddo
         atom_bas_done(iatom) = .true.

         ! Skips empty lines or those starting with #.
         line_read = ""
         do while (line_read == "")
            read(9999,*) line_read
            read(line_read, '(A1)') start_str
            if (start_str = "#") line_read = ""
         enddo
         read(line_read, '(A8)') start_str
      enddo
      close(9999)

      ! Reads fitting set data.
      ! Skips empty lines or those starting with #.
      open(unit = 9999, file= basis_file, iostat = file_iostat)
      line_read = ""
      do while (line_read == "")
         read(9999,*) line_read
         read(line_read, '(A1)') start_str
         if (start_str = "#") line_read = ""
      enddo
      read(line_read, '(A8)') start_str

      do while (start_str .ne. "endbasis")
         read(9999,*) iatom, nraw, ncon
         if (max_f_per_atom .lt. nraw) max_f_per_atom = nraw
         if (max_c_per_atom .lt. ncon) max_c_per_atom = ncon

         ! Only needs to read angular momenta for each contraction, the rest is
         ! skipped.
         do icount = 1, ncon
            read(9999,*)
         enddo
         do icount = 1, ncon
            read(9999,*) ang_mom
            if (any(atom_Z == iatom)) then
               aux_size = aux_size + ang_deg(ang_mom) * atom_count(iatom)
            endif
         enddo
         do icount = 1, nraw
            read(9999,*)
         enddo
         atom_fit_done(iatom) = .true.

         ! Skips empty lines or those starting with #.
         line_read = ""
         do while (line_read == "")
            read(9999,*) line_read
            read(line_read, '(A1)') start_str
            if (start_str = "#") line_read = ""
         enddo
         read(line_read, '(A8)') start_str
      enddo
      close(9999)
   endif
end subroutine basis_set_size

subroutine check_basis(atom_bas, atom_fit, n_atoms, atom_Z, iostatus)
   implicit none
   integer, intent(in)  :: n_atoms, atom_z(n_atoms)
   logical, intent(in)  :: atom_bas_done(0:120), atom_fit_done(0:120)
   integer, intent(out) :: iostat

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


! TAKE THIS TO OTHER MODULE
subroutine atom_name(atom_Z, symb)
 ! Takes atomic number Z and translates it to its name.
 implicit none
 integer         , intent(in)  :: atom_z
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
