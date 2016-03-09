######################################################################
# LIBLIO MAKEFILE
######################################################################
#
# DESIGN NOTES:
#
# The compilation process for the liblio-g2g.so library is currently
# being organized by this file and its companions: "make_depend.mk",
# "make_flags.mk" and "make_links.mk". The general scheeme is defined
# here, as well as all global/general variables to be used, and
# customization can be done by modifying the companion files. For
# the most part, this file should be left untouched (unless it is
# the compilation process itself what you want to modify).
#
# make_depend: This file contains the list of all objects and
#              modules to be compiled and incorporated in the
#              library. It also contains the information on how
#              these objects depend on one another and specifies
#              the paths to the sources, be those code files or .mk
#              files containing information on internal dependencies
#              for the modules on non compiling source files (which
#              are included in the compiling files instead).
#
# make_flags: This file contains the information on the flags to
#             be used during the compilaton of the objects. This
#             includes default flags as well as specific for a
#             group of objects or optional ones.
#
# make_links: This file contains the information on the flags to
#             be used during the final linking of the library, as
#             well as the information of the external libraries to
#             be linked (as well as other internal ones, such as
#             g2g which has to be precompiled).
#
######################################################################

# First things firs: we indicate the main target of the makefile.
all : liblio-g2g.so


# This variable will contain all .mk files that provide the details
# of how to compile the library. All objects will depend on this
# files, so that any modification in compilation options will
# trigger the recompilation of the code. We add here the basic ones
# and the rest should be incorporated to this variable in the
# make_depend file.
makefiles :=
makefiles += Makefile
makefiles += make_depend.mk
makefiles += make_flags.mk
makefiles += make_links.mk


# Variables are defined to list all the objects to be compiled and
# incorporated into liolib ("objects"), as well as the list of places
# where to find the source files ("src_path") and a single place
# where to put the objects and modules compiled ("obj_path").
# We set this paths for the makefile using vpath, which is not
# reliable for finding files that change during compilation time
# (that is, the .o and .mod) and thus it is not used for targets
# but for dependencies only, so the "obj_path" variable needs to
# be added explicitly in some places.
objects   :=
obj_path  := obj
src_paths := 
include make_depend.mk
vpath %.o   $(obj_path)
vpath %.mod $(obj_path)
vpath %.f   $(src_paths)
vpath %.f90 $(src_paths)
vpath %.mk  $(src_paths)

# This indicates which flags to use during the compilation
# of the different objects and which flags to use during
# the final linking of the library. It also contains the
# options to make with different compilers.
# LIBS is commented because for the time being LIO may rely
# on information set via environment variables.
FC =
FFLAGS =
LFLAGS =
#LIBS   =
include make_flags.mk
include make_links.mk


# GENERAL RULES
######################################################################
#   The following rules have been written to be general and make use
# of variables defined along this file and filled by the information
# provided by the other .mk files.
#   These rules should not be modified under normal circumstances;
# modifications of added files or change in the compilation flags
# and regular customization should be done through the variables
# defined before by modifying them inside the respective .mk file.

liblio-g2g.so : $(objects:%.o=$(obj_path)/%.o)
	$(FC) $(LFLAGS) $^ $(LIBS) -o $@

$(obj_path)/%.mod : $(obj_path)/%.o
	@touch $@

$(obj_path)/%.o   : %.f   $(makefiles) |  $(obj_path)
	$(FC) $(FFLAGS) -c $< -o $@

$(obj_path)/%.o   : %.F   $(makefiles) |  $(obj_path)
	$(FC) $(FFLAGS) -c $< -o $@

$(obj_path)/%.o   : %.f90 $(makefiles) |  $(obj_path)
	$(FC) $(FFLAGS) -c $< -o $@

$(obj_path)/%.o   : %.F90 $(makefiles) |  $(obj_path)
	$(FC) $(FFLAGS) -c $< -o $@

$(obj_path) :
	mkdir -p $@

.PHONY: clean
clean:
	rm -rf liblio-g2g.so $(obj_path)

######################################################################
