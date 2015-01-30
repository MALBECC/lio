######################################################################
# LIBLIO MAKEFILE
######################################################################

makefiles := Makefile
# This is the main file directing the compilation for the 
# library LIO. It was designed to simplify the customization
# of every aspect of this process and thus it should not be
# modified under normal code tinkering. Instead, one should
# attempt first to modify the auxiliary .mk files (beware
# not to change these files' names though).

makefiles += make_depend.mk
# This should contain the list of all objects and modules
# to be compiled and their internal dependencies/relations.
# It should also specify the locations of the source files
# (.f files) for these targets and include all .mk files
# containing information on how to compile modules.

makefiles += make_flags.mk
# This should contain the default flags for compilation
# as well as the specific flags for those targets that
# require them, and optional flags to be defined during
# compilation by the user.

makefiles += make_links.mk
# This should contain the default flags for final linking
# as well as the optional linking flags to be defined
# during compilation by the user. This also includes all
# external libraries that are to be linked to LIO.

#
######################################################################

all : liblio-g2g.so
# This indicates the main target.

obj_path := obj
vpath %.o   $(obj_path)
vpath %.mod $(obj_path)
# This indicates where to store the objects and modules produced
# by the compilation process. Vpath is not reliable for finding
# these files; obj_path shoud be added to objects explicitly.

objects   =
src_paths =
include make_depend.mk
vpath %.f  $(src_paths)
vpath %.mk $(src_paths)
# This indicates which objects are to be compiled, in which
# order, and where to look for its source files (or other
# dependencies, like .mk files)

FC     = ifort
#LIBS   =
FFLAGS =
LFLAGS =
include make_flags.mk
include make_links.mk
# This indicates which flags to use during the compilation
# of the different objects and which flags to use during
# the final linking of the library.
# LIBS is commented because for the time being LIO may rely
# on information set via environment variables.

#
######################################################################
# RULES
liblio-g2g.so: $(objects:%.o=$(obj_path)/%.o)
	$(FC) $(LFLAGS) $^ $(LIBS) -o $@

$(obj_path)/%.o   :  %.f $(makefiles) |  $(obj_path)
	$(FC) $(FFLAGS) -c $< -o $@ 

$(obj_path)/%.mod :  %.f $(makefiles) |  $(obj_path)
	$(FC) $(FFLAGS) -c $< -o $(@:%.mod=%.o)

$(obj_path) :
	mkdir -p $@

.PHONY: clean
clean:
	rm -rf liblio-g2g.so $(obj_path)

######################################################################
