!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  interface multiply_matrices
    module procedure matmult_rr
    module procedure matmult_dd
    module procedure matmult_cc
    module procedure matmult_zz

    module procedure matmult_rc
    module procedure matmult_cr
    module procedure matmult_dz
    module procedure matmult_zd

    module procedure matmult_rd
    module procedure matmult_dr
    module procedure matmult_cz
    module procedure matmult_zc

    module procedure matmult_rz
    module procedure matmult_zr
    module procedure matmult_dc
    module procedure matmult_cd
  end interface
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
