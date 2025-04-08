################################################################################
# Makefile - liblio
#
# The following makefile coordinates the compilation of the liblio library.
# This is a general purpose makefile, with generic rules that use generic
# variables. Please instead of modifying it directly, customize the process
# through the provided variables in the corresponding extra makefiles:
#
#
# (*) Makefile.depends contains information about which objects should be
#     compiled, how objects depend on one another, and where to look for the 
#     files to be compiled. Check this file when you are adding new objects.
#     This is where the variables OBJECTS and SRCDIRS are filled.
#
# (*) Makefile.options: contains information on the aplication of the
#     different compilation options. Check this file when you are customizing
#     how an object should be compiled.
#     This is where the variables FC, FFLAGS, LC, LIBS, LFLAGS are filled.
#
# (*) Makefile.optchck: script-like section to create and check a log on
#     how the code is being compiled. This allows to keep a record and keep
#     compilations consistent. Should not be modified frequently.
#
################################################################################

#.DEFAULT_GOAL := liblio-g2g.so
.DEFAULT_GOAL := all

FC =
FFLAGS =

LC =
LIBS =
LFLAGS =

OBJECTS :=
SRCDIRS :=
OBJPATH := obj
vpath %.o   $(OBJPATH)
vpath %.mod $(OBJPATH)

MAKEFILES =  Makefile
MAKEFILES += Makefile.optchck
MAKEFILES += Makefile.depends
MAKEFILES += Makefile.options

include Makefile.optchck
include Makefile.depends
include Makefile.options

vpath %.f   $(SRCDIRS)
vpath %.F   $(SRCDIRS)
vpath %.f90 $(SRCDIRS)
vpath %.F90 $(SRCDIRS)
vpath %.mk  $(SRCDIRS)

################################################################################
# RULES SECTION
#
# All objects depend on the important MAKEFILES because changes in the
# compilation process should always cause all objects to be recompiled.
# OBJPATH is set as an order-only prerequisite because we only want to
# make sure it is there, but we don't want to have every object recompile
# just because something was modified inside of it (else it would always
# cause full recompilation).

$(OBJPATH) :
	mkdir -p $@

%.mod : %.o             $(MAKEFILES) | $(OBJPATH)
	@touch $@

$(OBJPATH)/%.o :  %.f   $(MAKEFILES) | $(OBJPATH)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJPATH)/%.o :  %.f90 $(MAKEFILES) | $(OBJPATH)
	$(FC) $(FFLAGS) -c $< -o $@

LIBRARY=liblio-g2g.so
LIBRARY_LIBXC=compile_both
ORIG_LIBXC=$(libxc) # 
MAKE_FLAGS=$(filter-out -- w,$(MAKEFLAGS)) # Eliminate dashes from MAKEFLAGS
## Targets
.PHONY: all
all: $(if $(filter-out 0,$(ORIG_LIBXC)),$(LIBRARY_LIBXC),$(LIBRARY)) 

liblio-g2g.so : $(OBJECTS:%.o=$(OBJPATH)/%.o)
	@echo "Building liblio-g2g.so with g2g library with no libxc support."
	$(LC) $^ $(LFLAGS) $(LIBS_NO_LIBXC) $(LIBS) -o $@

liblio-g2g-libxc.so : $(OBJECTS:%.o=$(OBJPATH)/%.o)
	@echo "Building liblio-g2g-libxc.so with libg2g-libxc.so with libxc support."
	$(LC) $^ $(LFLAGS) $(LIBS_LIBXC) $(LIBS) -o $@

compile_both :
	@echo "Compiling with libxc=$(ORIG_LIBXC)..."
	$(MAKE) liblio-g2g-libxc.so libxc=$(ORIG_LIBXC) $(filter-out libxc=%, $(MAKE_FLAGS))
	@echo "Compiling with libxc=0..."
	$(MAKE) clean_obj
	$(MAKE) liblio-g2g.so libxc=0 $(filter-out libxc=%, $(MAKE_FLAGS))


.PHONY: .FORCE
.FORCE:

.PHONY: clean
clean:
	rm -rf Makefile.log liblio-g2g.so liblio-g2g-libxc.so *.o *.mod $(OBJPATH)

.PHONY: clean_obj
clean_obj:
	rm -rf Makefile.log *.o *.mod $(OBJPATH)

.PHONY: install
install:
	@echo 'This option is not available yet.'

.PHONY: uninstall
uninstall:
	@echo 'This option is not available yet.'

################################################################################
