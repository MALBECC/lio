!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module linear_algebra
!------------------------------------------------------------------------------!
  implicit none
  private
  integer, parameter :: SPK=selected_real_kind( 6, 37)
  integer, parameter :: DPK=selected_real_kind(15, 37)

  public :: Is_matrix_square
# include  "checks_interface.f90"

  public :: matrix_diagon
  interface matrix_diagon
    module procedure matrix_diagon_d
  end interface

  public :: multiply_matrices
# include "matmult_interface.f90"

  public :: commute_matrices
# include "matcommut_interface.f90"
!
!
!------------------------------------------------------------------------------!
contains
# include "is_matrix_square.f90"

# include "matrix_diagon_d.f90"
# include "matrix_diagon_dsyevd.f90"
# include "matrix_diagon_dsyevr.f90"

# include "matmult_procedures.f90"
# include "matcommut_procedures.f90"
end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
