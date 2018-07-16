.DEFAULT_GOAL := all

# -*- mode: make -*-
include Makefile.translate
include Makefile.common

## Define sources and objects
SRCS:=$(wildcard *.cpp)
OBJ:=$(SRCS:%.cpp=%.o)

POINTXC_SRC  := $(wildcard pointxc/*.cpp)
SRCS         += $(POINTXC_SRC)
OBJ          += $(POINTXC_SRC:%.cpp=%.o)

AINT_SRC     := $(wildcard analytic_integral/*.cpp)
SRCS         += $(AINT_SRC)
OBJ          += $(AINT_SRC:%.cpp=%.o)

CPU_SOURCES  := $(wildcard cpu/*.cpp)
SRCS         += $(CPU_SOURCES)
OBJ          += $(CPU_SOURCES:%.cpp=%.o)
CXXFLAGS     += -DCPU_KERNELS=1 -I$(MKLROOT)/include -Wno-unused-function
CPU          := 1

ifneq ($(cuda),0)
  include Makefile.cuda
  CUDA_LIBS    := -lcudart -lcuda
  SRCS         += $(CUDA_SOURCES)
  OBJ          += $(UGO) $(CUDA_OBJ)
  CXXFLAGS     += -DGPU_KERNELS=1
  NVCCFLAGS    += -DGPU_KERNELS=1
  LIBRARIES    := $(CUDA_LDFLAGS) $(CUDA_LIBS) -lrt
endif

ifeq ($(libxc),1)
  LIBXC_LIBRARIES := -L$(LIBXC_LIBS) -lxc
#LIBXC CPU MODE
  CXXFLAGS += -DUSE_LIBXC=1 -DLIBXC_CPU=1 -I$(LIBXC_INCLUDES)
else
ifeq ($(libxc),2)
#LIBXC GPU MODE
    LIBXC_LIBRARIES := -L$(LIBXC_LIBS) -lxc_cuda
    CXXFLAGS += -DUSE_LIBXC=1 -DLIBXC_CPU=0 -I$(LIBXC_INCLUDES)
else
#NO LIBXC
    LIBXC_LIBRARIES := -L$(LIBXC_LIBS) -lxc
    CXXFLAGS += -DUSE_LIBXC=0 -DLIBXC_CPU=1 -I$(LIBXC_INCLUDES)
endif
endif

ifneq ($(aint_gpu_level),)
  CXXFLAGS     += -DAINT_GPU_LEVEL=$(aint_gpu_level)
  NVCCFLAGS    += -DAINT_GPU_LEVEL=$(aint_gpu_level)
endif
ifeq ($(aint_mp),1)
  CXXFLAGS     += -DAINT_MP=1
  NVCCFLAGS    += -DAINT_MP=1
endif

## Other options
ifeq ($(print), 1)
        CXXFLAGS += -DPRINT_MATRICES
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
LIBRARIES := $(CUDA_LDFLAGS) $(CUDA_LIBS) $(LIBXC_LIBRARIES) -lrt
ifeq ($(intel),0)
  LIBRARIES += -L$(PACK_LIB_DIR) -lcblas
  CXXFLAGS += -DINTEL_COMP=0
else 
  CXXFLAGS += -DINTEL_COMP=1
endif


LIBRARY=libg2g.so
## Targets
all: $(LIBRARY)

include Makefile.depends

depend:
	@touch Makefile.depends
	@makedepend -Y -fMakefile.depends $(SRCS) > /dev/null 2>&1
	#@sed -ri 's|(cuda/.*)\.o|\1.cu_o|g' Makefile.depends
	@rm -f Makefile.depends.bak

libg2g.so: $(OBJ)
	$(CXX) -fopenmp -shared $(LDFLAGS) -o libg2g.so $(OBJ) $(LIBRARIES)

clean:
	@echo "Removing objects"; rm -f *.o libg2g.so *.a cpu/*.o cuda/*.o pointxc/*.o
	@rm -f cuda/*.cu_o cuda/*.cudafe* cuda/*.ptx cuda/*.hash cuda/*.cubin cuda/*.i cuda/*.ii cuda/*.fatbin.* cuda/*.cu.c
	@rm -f analytic_integral/*.o analytic_integral/cuda/*.o
