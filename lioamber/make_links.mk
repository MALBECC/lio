######################################################################
# LINK FLAGS : This are the flags related to the linkage of all
# the objects, modules and libraries (internal and external) to
# compile the LIO library.
LFLAGS += -shared -fPIC

#
######################################################################
# LIBRARIES : Included here are all the external libraries to be
# linked in the making of the LIO library (thogether with their
# respective paths). One can easily include new libraries by adding
# them (and their location) here.
LIBS += -L../g2g -L/usr/lib -L/usr/lib64
LIBS += -lg2g -lstdc++
LIBS += -liomp5 -lpthread -lm -Wl,-rpath='$$ORIGIN/' -Wl,-rpath='$$ORIGIN/../g2g'
LIBS += -L$(MKLROOT)/lib/intel64 -I$(MKLROOT)/include
LIBS += -lmkl_lapack95_lp64 -lmkl_intel_lp64
LIBS += -lmkl_intel_thread -lmkl_core

ifeq ($(magma),1)
  LIBS += -L$(MAGMAROOT)/lib -lmagma
endif

######################################################################
