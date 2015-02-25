######################################################################
# GENERAL FLAGS : This sets the configuration for the general
# flags that will be used in compiling all of the target objects
# (both for free subroutines and modules alike).
ifeq ($(ifort),1)
  FC     = ifort
  FFLAGS+= -module $(obj_path)
  FFLAGS+= -fpp
else
  FC     = gfortran
  FFLAGS+= -I$(obj_path) -J$(obj_path)
  FFLAGS+= -cpp
endif
FFLAGS += -Dpack -fPIC -DG2G -g
FFLAGS += $(optim) $(myflags) $(DEFINE) $(PROFILE)


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

objlist := dip.o SCFop.o
objlist += intfld.o int1G.o int3G.o
objlist += intsolG.o intsolGs.o intsol.o
$(objlist:%.o=$(obj_path)/%.o) : optim:=$(optim1)
#$(objlist:%.o=$(obj_path)/%.o) : private optim+=$(optim1)

objlist := matmuldiag.o int3lu.o
objlist := SCF.o TD.o ehrenfest.o magnus.o predictor.o
objlist += FixMessRho.o get_unit.o mulliken.o PackedStorage.f
objlist += init_amber.o init.o lio_init.o liomain.o lio_finalize.o
objlist += dft_get_mm_forces.o dft_get_qm_forces.o
objlist += alg.o drive.o func.o grid.o dipmem.o jarz.o
objlist += int1.o int2.o int2G.o int3mem.o int3mems.o intSG.o
objlist += garcha_mod.o mathsubs.o
$(objlist:%.o=$(obj_path)/%.o) : optim:=$(optim3)
#UNNECESSARY IF PREVIOUS ASSIGNMENT USED PRIVATE KEYWORD

#
######################################################################
# TARGET-SPECIFIC FLAGS (GENERAL) : This following section has
# all rest of the target-specific flags. More of these can be
# assigned to individual objects or to a whole list of them,
# grouped together inside the objlist, simultaneusly.
# (*)
myflags :=

ifeq ($(ifort),1)
  objlist := matmuldiag.o int3lu.o
  objlist += mathsubs.o
  $(objlist:%.o=$(obj_path)/%.o) : myflags:=-parallel
  #$(objlist:%.o=$(obj_path)/%.o) : private myflags+=-parallel

  objlist := SCF.o TD.o ehrenfest.o magnus.o predictor.o
  objlist += FixMessRho.o get_unit.o mulliken.o PackedStorage.f
  objlist += init_amber.o init.o lio_init.o liomain.o lio_finalize.o
  objlist += dft_get_mm_forces.o dft_get_qm_forces.o
  objlist += alg.o drive.o func.o grid.o dipmem.o jarz.o
  objlist += int1.o int2.o int2G.o int3mem.o int3mems.o intSG.o
  objlist += garcha_mod.o
  $(objlist:%.o=$(obj_path)/%.o) : myflags:=-mp1 -ip
  #$(objlist:%.o=$(obj_path)/%.o) : private myflags+=$(optim3) -mp1 -ip
endif
#
######################################################################
# (*) IMPORTANT: due to the incompatibility of most common
# versions of make with the 'private' keyword, the way target-
# -specific flags are implemented here can lead to problems
# related to the inheritance of these flags; for example, the
# lack of reliability in any default value for optimization.
######################################################################
