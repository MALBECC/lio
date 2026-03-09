#include "../datatypes/datatypes.fh"
module density_fitting_verbosity
  implicit none
  private
  public :: write_af_record, write_propd

contains

  ! Return a free I/O unit number (F90-safe)
  integer function get_free_unit()
    implicit none
    integer :: u
    logical :: is_open
    do u = 10, 999
      inquire(unit=u, opened=is_open)
      if (.not. is_open) then
        get_free_unit = u
        return
      end if
    end do
    get_free_unit = -1
  end function get_free_unit


  ! Write AF block
  subroutine write_af_record(niter, Md, E1, E2, kinE, Exc, En0, af, &
                             filename, append_mode)
    implicit none
    ! inputs
    integer, intent(in) :: niter, Md
    LIODBLE, intent(in) :: E1, E2, kinE, Exc, En0
    LIODBLE, intent(in) :: af(Md)
    character(*), intent(in), optional :: filename
    logical, intent(in),     optional :: append_mode
    ! locals
    integer :: u, ios, i
    character(len=256) :: fname
    logical :: do_append

    fname = 'df_coeffs.log'
    if (present(filename)) fname = filename
    do_append = .true.
    if (present(append_mode)) do_append = append_mode

    u = get_free_unit()
    if (u .lt. 0) then
      write(*,*) 'density_fitting_verbosity: Error: no free unit available to open "', trim(fname), '"'
      return
    end if

    if (do_append) then
      open(unit=u, file=trim(fname), status='unknown', position='append', &
           action='write', iostat=ios)
    else
      open(unit=u, file=trim(fname), status='replace', action='write', iostat=ios)
    end if
    if (ios .ne. 0) then
      write(*,*) 'density_fitting_verbosity: Error opening file "', trim(fname), '", iostat=', ios
      return
    end if

    write(u, *) niter
    write(u, *) Md
    write(u, *) E1
    write(u, *) E2
    write(u, *) kinE
    write(u, *) Exc
    write(u, *) En0 + E1 + E2 + Exc

    do i = 1, Md
      write(u, '(F18.9)') af(i)
    end do
    write(u, *) 'END'

    close(u)
  end subroutine write_af_record


  ! Write PROPD (nuc index + first components of ad, cd + angular momentum + contractions)
  subroutine write_propd(Nucd, ad, cd, Md, ang_momd, nContd, filename, append_mode)
    implicit none
    integer, intent(in) :: Md
    integer, intent(in) :: Nucd(Md)
    integer, intent(in) :: ang_momd(Md)
    integer, intent(in) :: nContd(Md)
    LIODBLE, intent(in) :: ad(Md, *), cd(Md, *)
    character(*), intent(in), optional :: filename
    logical, intent(in),     optional :: append_mode
    integer :: u, ios, i
    character(len=256) :: fname
    logical :: do_append

    fname = 'df_props.log'
    if (present(filename)) fname = filename
    do_append = .true.
    if (present(append_mode)) do_append = append_mode

    u = get_free_unit()
    if (u .lt. 0) then
      write(*,*) 'density_fitting_verbosity: Error: no free unit available to open "', trim(fname), '"'
      return
    end if

    if (do_append) then
      open(unit=u, file=trim(fname), status='unknown', position='append', &
           action='write', iostat=ios)
    else
      open(unit=u, file=trim(fname), status='replace', action='write', iostat=ios)
    end if
    if (ios .ne. 0) then
      write(*,*) 'density_fitting_verbosity: Error opening file "', trim(fname), '", iostat=', ios
      return
    end if

    do i = 1, Md
      write(u,*) Nucd(i), ang_momd(i), nContd(i), ad(i,1), cd(i,1)
    end do

    close(u)
  end subroutine write_propd

end module density_fitting_verbosity

