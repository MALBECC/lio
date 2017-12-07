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

.DEFAULT_GOAL := liblio-g2g.so

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

liblio-g2g.so : $(OBJECTS:%.o=$(OBJPATH)/%.o)
	$(LC) $^ $(LFLAGS) $(LIBS) -o $@

.PHONY: .FORCE
.FORCE:

.PHONY: clean
clean:
	rm -rf Makefile.log liblio-g2g.so *.o *.mod $(OBJPATH)

.PHONY: install
install:
	@echo 'This option is not available yet.'

.PHONY: uninstall
uninstall:
	@echo 'This option is not available yet.'

################################################################################
