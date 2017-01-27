######################################################################
# GENERAL FLAGS : This sets the configuration for the general
# flags that will be used in compiling all of the target objects
# (both for free subroutines and modules alike).
ifeq ($(intel),1)
  FC     = ifort
  FFLAGS+= -module $(obj_path)
  FFLAGS+= -fpp
#  FFLAGS+= -check bounds
  FFLAGS+= -traceback -check all -fp-stack-check
else
  FC     = gfortran
  FFLAGS+= -I$(obj_path) -J$(obj_path)
  FFLAGS+= -cpp
endif
DEFINE  += -Dpack -DG2G
PROFILE += -g
FFLAGS  += -fPIC
FFLAGS  += $(DEFINE) $(PROFILE) $(optim) $(myflags)


#
######################################################################
# OPTIONAL FLAGS : This sets the different optional flags to
# be set via the command line when calling  make to compile the
# code. Note that the options on how to optimize compilation are
# dealt with in the next section.
ifeq ($(print), 1)
  DEFINE += -DPRINT_MATRICES
endif

ifeq ($(profile),1)
  PROFILE = -pg
else
  PROFILE =
endif

ifeq ($(magma),1)
  DEFINE += -Dmagma
endif

FORTRAN_WRAPPER=
ifeq ($(cublas),1)
  DEFINE += -DCUBLAS
endif

ifeq ($(td_simple),1)
  DEFINE += -DTD_SIMPLE
endif


# Only uncomment this if your version of Amber supports
# lio/Amber basis set keywords and cube file generation
#FFLAGS += -DMOD_AMBER

#
######################################################################
# TARGET-SPECIFIC FLAGS (OPTIMIZATION) : This following section
# assigns to each object to be compiled the optimization with
# which to do it.
# (*)
optim1=-O1
optim3=-O3
ifeq ($(non_optimize),1)
  optim1=-O0
  optim3=-O0
else ifeq ($(full_optimize),1)
  optim1=-O3
  optim3=-O3
endif
optim=$(optim3)

tmplist := dip.o SCFop.o
tmplist += intfld.o int1G.o int3G.o
tmplist += intsolG.o intsolGs.o intsol.o
$(tmplist:%.o=$(obj_path)/%.o) : optim:=$(optim1)
#$(tmplist:%.o=$(obj_path)/%.o) : private optim+=$(optim1)

tmplist := matmuldiag.o int3lu.o fock_commuts.o
tmplist := SCF.o TD.o ehrenfest.o magnus.o predictor.o
tmplist += FixMessRho.o mulliken.o PackedStorage.f
tmplist += init_amber.o init.o lio_init.o liomain.o lio_finalize.o
tmplist += dft_get_mm_forces.o dft_get_qm_forces.o
tmplist += alg.o drive.o func.o grid.o dipmem.o jarz.o
tmplist += int1.o int2.o int2G.o int3mem.o intSG.o
tmplist += garcha_mod.o mathsubs.o cubegen.o density.o

ifeq ($(cublas),1)
tmplist += cublasmath.o 
endif
$(tmplist:%.o=$(obj_path)/%.o) : optim:=$(optim3)
#UNNECESSARY IF PREVIOUS ASSIGNMENT USED PRIVATE KEYWORD

#
######################################################################
# TARGET-SPECIFIC FLAGS (GENERAL) : This following section has
# all rest of the target-specific flags. More of these can be
# assigned to individual objects or to a whole list of them,
# grouped together inside the tmplist, simultaneusly.
# (*)
myflags :=

ifeq ($(intel),1)
  tmplist := matmuldiag.o int3lu.o fock_commuts.o
  tmplist += mathsubs.o
  $(tmplist:%.o=$(obj_path)/%.o) : myflags:=-parallel
  #$(tmplist:%.o=$(obj_path)/%.o) : private myflags+=-parallel

  tmplist := SCF.o TD.o ehrenfest.o magnus.o predictor.o
  tmplist += FixMessRho.o mulliken.o PackedStorage.f
  tmplist += init_amber.o init.o lio_init.o liomain.o lio_finalize.o
  tmplist += dft_get_mm_forces.o dft_get_qm_forces.o
  tmplist += alg.o drive.o func.o grid.o dipmem.o jarz.o
  tmplist += int1.o int2.o int2G.o int3mem.o  intSG.o
  tmplist += garcha_mod.o cubegen.o density.o
  $(tmplist:%.o=$(obj_path)/%.o) : myflags:=-mp1 -ip
  #$(tmplist:%.o=$(obj_path)/%.o) : private myflags+=$(optim3) -mp1 -ip
endif


#
######################################################################
# (*) IMPORTANT: due to the incompatibility of most common
# versions of make with the 'private' keyword, the way target-
# -specific flags are implemented here can lead to problems
# related to the inheritance of these flags; for example, the
# lack of reliability in any default value for optimization.
######################################################################
