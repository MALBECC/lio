! This file is part of the LIO software package (https://github.com/lio-project/lio)
! This is the part of liosolo with support for libxc functionality. 
! This is file should be linked against liblio-g2g.so.
subroutine liosolo_native(escf, dipxyz)
    use garcha_mod , only : natom, nsol, Iz, writeforces, writexyz
    use basis_data , only : basis_set, fitting_set, int_basis
    use ECP_mod    , only : ecpmode
    use fileio_data, only : verbose
    use fileio     , only : lio_logo

    implicit none
    real*8            :: dipxyz(3), escf

    ! Initialize lio with libxc. 
    call init_lio_common(natom, Iz, nsol, 0)

    ! Calls main procedures.
    call liomain(escf, dipxyz)
    call lio_finalize()
end subroutine liosolo_native
