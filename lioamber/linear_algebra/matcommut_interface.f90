!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  interface commute_matrices
    module procedure matcommut_rr
    module procedure matcommut_dd
    module procedure matcommut_cc
    module procedure matcommut_zz

    module procedure matcommut_rc
    module procedure matcommut_cr
    module procedure matcommut_dz
    module procedure matcommut_zd

    module procedure matcommut_rd
    module procedure matcommut_dr
    module procedure matcommut_cz
    module procedure matcommut_zc

    module procedure matcommut_rz
    module procedure matcommut_zr
    module procedure matcommut_dc
    module procedure matcommut_cd
  end interface
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
