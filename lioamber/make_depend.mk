######################################################################
# SOURCE LOCATIONS : This is the list of all the places where
# to look for the source files to compile the targets.
src_paths += liomods
src_paths += maskrmm
src_paths += mathsubs

#
######################################################################
# MODULE INTERNALS : This is the inclusion of all .mk files
# that contain the internal information on how to compile the
# different modules.
include liomods/liomods.mk
include maskrmm/maskrmm.mk
include mathsubs/mathsubs.mk

#
######################################################################
# OBJECTS LIST : This is the list of ALL the objects to be
# linked for making LIO, including both free subroutines and
# modules (their corresponding .o files).
objects += liomain.o SCF.o SCFop.o SCF_in.o TD.o cubegen.o
objects += dip.o dipmem.o jarz.o magnus.o predictor.o mulliken.o
objects += dft_get_mm_forces.o dft_get_qm_forces.o
objects += matmuldiag.o fock_commuts.o
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
objects += elec.o

#
######################################################################
# OTHER SPECIFICS : This is the list of all the specific
# dependencies, including dependence on modules and on other
# separated files. These can be specified individually or by
# initializing the objlist. Please note that the obj_path is
# explicitly set in each case; vpath can't be trusted with
# handling files generated during compiling.


# garcha_mod
objlist := dft_get_mm_forces.o dft_get_qm_forces.o
objlist += dip.o dipmem.o drive.o grid.o init_amber.o init.o
objlist += int1.o int2.o int3lu.o int3mem.o int3mems.o intfld.o intsol.o
objlist += int1G.o int2G.o int3G.o intSG.o intsolG.o intsolGs.o
objlist += jarz.o lio_finalize.o predictor.o
objlist += SCF.o SCF_in.o SCFop.o TD.o cubegen.o
$(objlist:%.o=$(obj_path)/%.o) : $(obj_path)/garcha_mod.o
$(objlist:%.o=$(obj_path)/%.o) : $(obj_path)/garcha_mod.mod


# mathsubs
objlist := SCF.o SCFop.o
$(objlist:%.o=$(obj_path)/%.o) : $(obj_path)/mathsubs.o
$(objlist:%.o=$(obj_path)/%.o) : $(obj_path)/mathsubs.mod

######################################################################
