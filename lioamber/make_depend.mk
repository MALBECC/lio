######################################################################
# This file has the dependency scheme of all lio files. It is
# the most likely to be modified, since any new objects to be
# included must be specified here.
#
# The file has "three" sections:
#
# (1) First a section to list all objects to be made, including
#     the '.o' files asociated to modules.
#
# (2) Then a section to specify the dependency on modules for
#     the objects of lio.
#
# (3) Finally a section to set specific flags for some of the
#     objects to be compiled.
#
# The way the first section works should be self evident.
# The way the second and third section work is by first starting
# a list of all objects to which apply the dependency or flag
# ("objlist" is the generic variable used for this purpose; it
# should be re-started again before usage) and then applying
# the modification to the whole list at once (and replacing
# with the right path; note that this variable is set in the
# main makefile).
#
#
######################################################################
# IMPORTANT: the target specific flags section has been
# implemented in a way in which inheritance of these variables
# is allowed, which has two consequences: (1) groups with no
# specification in them may inherit one because of being the
# dependency of another target that has it; (2) in order to
# override inherited variables, "myflag" has to be overwritten
# instead of just appending the flags, so all assignment groups
# have to be excluding (targets should not appear in two groups).
#
# This current inconvinience is due to the incompatibility of the
# keyword 'private' with the version of GNUmake included in the
# most common repositories. This should be open to review in the
# near future (either because GNUmake version 3.82 has become
# commonplace or a more flexible non-exclusive grouping method
# has become more convenient for target-specific assignment).
#
#
######################################################################
# Object List
objects += liomain.o SCF.o SCFop.o TD.o
objects += dip.o dipmem.o jarz.o magnus.o predictor.o mulliken.o
objects += dft_get_mm_forces.o dft_get_qm_forces.o
objects += matmuldiag.o
objects += init.o init_amber.o lio_init.o lio_finalize.o
objects += alg.o drive.o func.o grid.o
objects += int1.o int1G.o int2.o int2G.o
objects += int3lu.o int3mem.o int3mems.o int3G.o
objects += intsol.o intsolG.o intsolGs.o
objects += intfld.o intSG.o
objects += FixMessRho.o get_unit.o PackedStorage.f
objects += garcha_mod.o
objects += liokeys.o
objects += sysdata.o
objects += mathsubs.o
objects += maskrmm.o
#
######################################################################
# garcha_mod
objlist := dft_get_mm_forces.o dft_get_qm_forces.o
objlist += dip.o dipmem.o drive.o grid.o init_amber.o init.o
objlist += int1.o int2.o int3lu.o int3mem.o int3mems.o intfld.o intsol.o
objlist += int1G.o int2G.o int3G.o intSG.o intsolG.o intsolGs.o
objlist += jarz.o lio_finalize.o predictor.o SCF.o SCF_in.o SCFop.o

$(objlist:%.o=$(obj_path)/%.o) : $(obj_path)/garcha_mod.o
$(objlist:%.o=$(obj_path)/%.o) : $(obj_path)/garcha_mod.mod
######################################################################
# mathsubs
objlist := SCF.o SCFop.o

$(objlist:%.o=$(obj_path)/%.o) : $(obj_path)/mathsubs.o
$(objlist:%.o=$(obj_path)/%.o) : $(obj_path)/mathsubs.mod
######################################################################
# Target specific flags
myflags :=
ifeq ($(non_optimize),1)
  optim1=-O0
  optim3=-O0
else
  optim1=-O1
  optim3=-O3
endif

objlist := matmuldiag.o int3lu.o
objlist += mathsubs.o
$(objlist:%.o=$(obj_path)/%.o) : myflags:=$(optim3) -parallel
#$(objlist:%.o=$(obj_path)/%.o) : private myflags+=$(optim3) -parallel

objlist := dip.o SCFop.o
objlist += intfld.o int1G.o int3G.o intsolG.o intsolGs.o intsol.o
$(objlist:%.o=$(obj_path)/%.o) : myflags:=$(optim1)
#$(objlist:%.o=$(obj_path)/%.o) : private myflags+=$(optim1)

objlist := SCF.o TD.o ehrenfest.o magnus.o predictor.o
objlist += FixMessRho.o get_unit.o mulliken.o PackedStorage.f
objlist += init_amber.o init.o lio_init.o liomain.o lio_finalize.o
objlist += dft_get_mm_forces.o dft_get_qm_forces.o
objlist += alg.o drive.o func.o grid.o dipmem.o jarz.o
objlist += int1.o int2.o int2G.o int3mem.o int3mems.o intSG.o
objlist += garcha_mod.o
$(objlist:%.o=$(obj_path)/%.o) : myflags:=$(optim3) -mp1 -ip
#$(objlist:%.o=$(obj_path)/%.o) : private myflags+=$(optim3) -mp1 -ip
######################################################################
