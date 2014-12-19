######################################################################
# LIBLIO MAKEFILE
######################################################################
# Main file directing the compilation process of lio. This
# file needs only be modified if a new folder with files is
# added, in which case one might need to:
#
# (1) Add the new folder to the src_paths so make checks this
#     location when searching for source .f files.
#
# (2) If there are internal dependencies of files that are
#     "hidden"*, these may be indicated in a ".mk" file and
#      then included here. This is how we can currently deal
#      with big modules.
#
# * (some files that may not compile directly but are included
#   in other files, and hence modifications on them would affect
#   the including file) 
#
#
# Other modifications should be more unlikely, and in case of
# necessity, should be done with care.
#
######################################################################
obj_path  := obj
src_paths := liomods
src_paths += maskrmm
src_paths += mathsubs
vpath %.mod $(obj_path)
vpath %.o   $(obj_path)
vpath %.f   $(src_paths)

include make_options.mk

all: liblio-g2g.so

objects :=
include make_depend.mk
include maskrmm/maskrmm.mk
include mathsubs/mathsubs.mk

liblio-g2g.so: $(objects:%.o=$(obj_path)/%.o)
	$(FC) -shared -fPIC -o $@ $^ $(LIBS)

include make_rules.mk
$(objects:%.o=$(obj_path)/%.o) : Makefile
$(objects:%.o=$(obj_path)/%.o) : make_rules.mk
$(objects:%.o=$(obj_path)/%.o) : make_depend.mk
$(objects:%.o=$(obj_path)/%.o) : make_options.mk
######################################################################
