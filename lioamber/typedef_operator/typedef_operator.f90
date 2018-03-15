module typedef_operator

   type operator
     real*8, allocatable :: data_AO(:,:)
     real*8, allocatable :: data_ON(:,:)

   contains
     procedure, pass :: Sets_data_AO
     procedure, pass :: Gets_data_AO
     procedure, pass :: Sets_data_ON
     procedure, pass :: Gets_data_ON
     procedure, pass :: Diagon_datamat
     procedure, pass :: Dens_build
     procedure, pass :: Commut_data
     procedure, pass :: BChange_AOtoON
   end type operator

contains
#include  "Sets_datamat.f90"
#include  "Gets_datamat.f90"
#include  "Diagon_datamat.f90"
#include  "Dens_build.f90"
#include  "Commut_data.f90"
#include  "BChange_data.f90"

end module typedef_operator
