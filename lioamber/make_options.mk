######################################################################
# This file contains the configuration of important variables
# in the compiling process of lio so it should be modified
# with due care.
#
# FLAGS: The flags here defined will be used for the compilation
# of all objects, so one should be careful when modifying it. If 
# some files require some specific flags, those can be defined
# privately in make_depend using the myflags variable.
#
# LIBS: Included here are the libraries and the respective paths
# where they reside. This variable is used in final linking with
# liolibs. One can easily include new libraries by adding them
# (and their location) here.
#
# OPTIONS: Options to include via command line specifications
# go here in the form of if conditionals. We recomend that this
# options be restricted to adding libraries or flags; other sort
# of customizations may require for upkeep to check modifications
# in the other make-related files (.mk) and is not advisable. Should
# such customizations be needed, a more drastic change in the make
# scheme would probably be advisable.
#
######################################################################
# FLAGS FOR ALL OBJECTS
FC     = ifort
FFLAGS =
FFLAGS += -Dpack -fPIC -DG2G
FFLAGS += -g -fpp -module $(obj_path)
FFLAGS += $(myflags) $(DEFINE) $(PROFILE)
######################################################################
# EXTERNAL LIBS
LIBS += -L../g2g -L/usr/lib -L/usr/lib64
LIBS += -lg2g -lstdc++
LIBS += -liomp5 -lpthread -lm -Wl,-rpath='$$ORIGIN/' -Wl,-rpath='$$ORIGIN/../g2g'
LIBS += -L$(MKLROOT)/lib/intel64 -I$(MKLROOT)/include
LIBS += -lmkl_lapack95_lp64 -lmkl_intel_lp64
LIBS += -lmkl_intel_thread -lmkl_core
######################################################################
# COMPILATION OPTIONS
ifeq ($(print), 1)
  DEFINE += -DPRINT_MATRICES
endif
#
ifeq ($(profile),1)
  PROFILE = -pg
else
  PROFILE =
endif
#
ifeq ($(magma),1)
  DEFINE += -Dmagma
  LIBS   += -L$(MAGMAROOT)/lib -lmagma
endif
######################################################################
