# -*- mode: make -*-
include Makefile.common

## Define sources and objects
SRCS:=$(wildcard *.cpp)
OBJ:=$(SRCS:%.cpp=%.o)

CPU_SOURCES  := $(wildcard cpu/*.cpp)
SRCS         += $(CPU_SOURCES)
OBJ          += $(CPU_SOURCES:%.cpp=%.o)
CXXFLAGS     += -DCPU_KERNELS=1 -I$(MKLROOT)/include
CPU          := 1

ifneq ($(cuda),0)
  include Makefile.cuda
  CUDA_LIBS    := -lcudart -lcuda
  CUDA_SOURCES := $(wildcard cuda/*.cu)
  SRCS         += $(CUDA_SOURCES)
  OBJ          += $(CUDA_SOURCES:%.cu=%.cu_o)
  CXXFLAGS     += -DGPU_KERNELS=1
  NVCCFLAGS    += -DGPU_KERNELS=1
  LIBRARIES    := $(CUDA_LDFLAGS) $(CUDA_LIBS) -lrt
endif

ifeq ($(coulomb_gpu),1)
  CXXFLAGS     += -DGPU_COULOMB=1
endif

## Other options
ifeq ($(print), 1)
        CXXFLAGS += -DPRINT_MATRICES
endif

ifeq ($(time), 1)
  CXXFLAGS += -DTIMINGS
endif

ifeq ($(profile),1)
  CXXFLAGS += -pg -g
  LDFLAGS += -pg -g
endif

ifeq ($(cpu_recompute),1)
  CXXFLAGS += -DCPU_RECOMPUTE=1
else
  CXXFLAGS += -DCPU_RECOMPUTE=0
endif

ifeq ($(full_double),1)
  CXXFLAGS += -DFULL_DOUBLE=1
  NVCCFLAGS += -DFULL_DOUBLE=1
else
  CXXFLAGS += -DFULL_DOUBLE=0
  NVCCFLAGS += -DFULL_DOUBLE=0
endif

## Define libraries
LIBRARIES := $(CUDA_LDFLAGS) $(CUDA_LIBS) -lrt
ifneq ($(ifort),1)
  LIBRARIES += -L$(PACK_LIB_DIR) -lcblas
endif

LIBRARY=libg2g.so
## Targets
all: $(LIBRARY)

include Makefile.depends

depend:
	@touch Makefile.depends
	@makedepend -Y -fMakefile.depends $(SRCS) > /dev/null 2>&1
	@sed -ri 's|(cuda/.*)\.o|\1.cu_o|g' Makefile.depends
	@rm -f Makefile.depends.bak

libg2g.so: $(OBJ)
	#@if [ -z "$(CUDA_HOME)" ] ; then echo "CUDA_HOME is not set. Please, check it."; exit 1; fi;
	#@if ! test -f $(CUDA_HOME)/lib/libcudart.so; then echo "libcudart.so can not be found in CUDA_HOME/lib/, current value is '"$(CUDA_HOME)/lib"'. Please check if CUDA_HOME variable is correctly set."; exit 1; fi;
	$(CXX) -fopenmp -shared $(LDFLAGS) -o libg2g.so $(OBJ) $(LIBRARIES)

clean:
	@echo "Removing objects"; rm -f *.o libg2g.so *.a cpu/*.o cuda/*.o
	@rm -f cuda/*.cu_o cuda/*.cudafe* cuda/*.ptx cuda/*.hash cuda/*.cubin cuda/*.i cuda/*.ii cuda/*.fatbin.* cuda/*.cu.c
