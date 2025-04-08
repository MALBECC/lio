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

# Define LIBXC CPU/GPU
LIBXC_CPU_LIBS=$(LIBXC_HOME_CPU)/lib
LIBXC_GPU_LIBS=$(LIBXC_HOME_GPU)/lib
LIBXC_CPU_INCLUDES=$(LIBXC_HOME_CPU)/include
LIBXC_GPU_INCLUDES=$(LIBXC_HOME_GPU)/include

ifneq ($(cuda),0)
  include Makefile.cuda
  CUDA_LIBS    := -lcudart -lcuda
  SRCS         += $(CUDA_SOURCES)
  OBJ          += $(UGO) $(CUDA_OBJ)
  CXXFLAGS     += -DGPU_KERNELS=1
  NVCCFLAGS    += -DGPU_KERNELS=1
  LIBRARIES    := $(CUDA_LDFLAGS) $(CUDA_LIBS) -lrt
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
##  LIBRARIES += -L$(PACK_LIB_DIR) -lcblas
  CXXFLAGS += -DINTEL_COMP=0
else 
  CXXFLAGS += -DINTEL_COMP=1
endif

ifeq ($(libxc),1)
  LIBXC_LIBRARIES := -L$(LIBXC_CPU_LIBS) -lxc
#LIBXC CPU MODE
  CXXFLAGS += -DUSE_LIBXC=1 -DLIBXC_CPU=1 -I$(LIBXC_CPU_INCLUDES)
else
ifeq ($(libxc),2)
#LIBXC GPU MODE
    LIBXC_LIBRARIES := -L$(LIBXC_GPU_LIBS) -lxc_cuda -L$(LIBXC_CPU_LIBS) -lxc
    CXXFLAGS += -DUSE_LIBXC=1 -DLIBXC_CPU=0 -I$(LIBXC_GPU_INCLUDES) -I$(LIBXC_CPU_INCLUDES)
else
#NO LIBXC
    LIBXC_LIBRARIES := -L$(LIBXC_LIBS) -lxc
    CXXFLAGS += -DUSE_LIBXC=0 -DLIBXC_CPU=1 -I$(LIBXC_INCLUDES)
endif
endif


include Makefile.depends
ifeq ($(libint),1)
    LIBS_LIBINT = -L/$(LIBINT_HOME)/lib -lint2
    INCS_LIBINT = -I/$(LIBINT_HOME)/include -I/$(LIBINT_HOME)/include/libint2 -I/$(EIGEN_HOME)
    LIBRARIES  += $(LIBS_LIBINT)
    CXXFLAGS   += -std=c++11 $(INCS_LIBINT) -DUSE_LIBINT=1

    SRCS       += libint/libintproxy.cpp libint/g2g_libint.cpp
    OBJ        += libint/libintproxy.o libint/g2g_libint.o
else
    SRCS       += libint/g2g_libint_fake.cpp
    OBJ        += libint/g2g_libint_fake.o
    CXXFLAGS   += -DUSE_LIBINT=0
endif
include Makefile.excited
include Makefile.cioverlap

LIBRARY=libg2g.so
LIBRARY_LIBXC=compile_both
ORIG_LIBXC=$(libxc) # 
MAKE_FLAGS=$(filter-out -- w,$(MAKEFLAGS)) # Eliminate dashes from MAKEFLAGS

## Targets
all: $(if $(filter-out 0,$(ORIG_LIBXC)),$(LIBRARY_LIBXC),$(LIBRARY)) 

depend:
	@touch Makefile.depends
	@makedepend -Y -fMakefile.depends $(SRCS) > /dev/null 2>&1
	#@sed -ri 's|(cuda/.*)\.o|\1.cu_o|g' Makefile.depends
	@rm -f Makefile.depends.bak

## Make library with libxc=1 or libxc=2
libg2g_libxc.so: $(OBJ) #libg2g.so
	@echo "Building libg2g_libxc.so library with LIBXC=1 or 2."
	$(CXX) -fopenmp -shared $(LDFLAGS) -o libg2g_libxc.so $(OBJ) $(LIBRARIES)

## Make library with libxc=0
libg2g.so: $(OBJ)
	@echo "Building libg2g.so library with LIBXC=0."
	$(CXX) -fopenmp -shared $(LDFLAGS) -o libg2g.so $(OBJ) $(LIBRARIES)

compile_both:
	@echo "Compiling first with libxc=$(ORIG_LIBXC)..."
	@echo "$(MAKE_FLAGS)"
	$(MAKE) libg2g_libxc.so libxc=$(ORIG_LIBXC) $(filter-out libxc=%, $(MAKE_FLAGS))
	@echo "Compiling second with libxc=0..."
	@echo "$(MAKE_FLAGS)"
	$(MAKE) clean_obj
	$(MAKE) libg2g.so libxc=0 $(filter-out libxc=%, $(MAKE_FLAGS))

clean:
	@echo "Removing objects"; rm -f *.o* libg2g.so libg2g_libxc.so *.a cpu/*.o* cuda/*.o* pointxc/*.o*
	@rm -f cuda/*.cu_o cuda/*.cudafe* cuda/*.ptx cuda/*.hash cuda/*.cubin cuda/*.i cuda/*.ii cuda/*.fatbin.* cuda/*.cu.c
	@rm -f analytic_integral/*.o* analytic_integral/cuda/*.o*
	@rm -f excited/*.o* excited/cuda/*.o* libint/*o

clean_obj:
	@echo "Cleaning objects for libxc build"
	rm -f *.o* cpu/*.o* cuda/*.o* pointxc/*.o*
	rm -f cuda/*.cu_o cuda/*.cudafe* cuda/*.ptx cuda/*.hash cuda/*.cubin cuda/*.i cuda/*.ii cuda/*.fatbin.* cuda/*.cu.c
	rm -f analytic_integral/*.o* analytic_integral/cuda/*.o*
	rm -f excited/*.o* excited/cuda/*.o* libint/*o
