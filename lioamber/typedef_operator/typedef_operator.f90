module typedef_operator
#include "../complex_type.fh"
   type operator
     real*8, allocatable     :: data_AO(:,:)
     real*8, allocatable     :: data_ON(:,:)
     TDCOMPLEX, allocatable  :: dataC_AO(:,:)
     TDCOMPLEX, allocatable  :: dataC_ON(:,:)

   contains
     procedure, pass :: Sets_data_AO
     procedure, pass :: Gets_data_AO
     procedure, pass :: Sets_dataC_AO
     procedure, pass :: Gets_dataC_AO
     procedure, pass :: Sets_data_ON
     procedure, pass :: Gets_data_ON
     procedure, pass :: Sets_dataC_ON
     procedure, pass :: Gets_dataC_ON
     procedure, pass :: Diagon_datamat
     procedure, pass :: Dens_build
     procedure, pass :: Commut_data_r
     procedure, pass :: Commut_data_c
     procedure, pass :: BChange_AOtoON
     procedure, pass :: BChange_ONtoAO
     procedure, pass :: Shift_diag_ON
     procedure, pass :: Shift_diag_AO
   end type operator

contains
#include  "Sets_datamat.f90"
#include  "Gets_datamat.f90"
#include  "Diagon_datamat.f90"
#include  "Dens_build.f90"
#include  "Commut_data.f90"
#include  "BChange_data.f90"
#include  "Shift_diag.f90"

end module typedef_operator
