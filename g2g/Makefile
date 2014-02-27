# -*- mode: make -*-
include Makefile.common

## Define sources and objects
SRCS:=$(wildcard *.cpp)
OBJ:=$(SRCS:%.cpp=%.o)

ifeq ($(cpu),1)
  include Makefile.cuda
  CPU_SOURCES  := $(wildcard cpu/*.cpp)
  SRCS         += $(CPU_SOURCES)
  OBJ          += $(CPU_SOURCES:%.cpp=%.o)
  CXXFLAGS     += -DCPU_KERNELS=1 -L/usr/local/cuda/cuda/lib64
else
  include Makefile.cuda
  CUDA_LIBS    :=-L/usr/local/cuda/cuda/lib64 -lcudart
  CUDA_SOURCES := $(wildcard cuda/*.cu)
  SRCS         += $(CUDA_SOURCES)
  OBJ          += $(CUDA_SOURCES:%.cu=%.cu_o)
endif

## Other options
ifeq ($(time), 1)
  CXXFLAGS += -DTIMINGS
endif

ifeq ($(nosync), 1)
  CXXFLAGS += -DNOSYNC
endif

ifeq ($(histo), 1)
  CXXFLAGS += -DHISTOGRAM
endif

ifeq ($(profile),1)
  CXXFLAGS += -pg
  LDFLAGS += -pg
endif

ifeq ($(cpu_recompute),0)
  CXXFLAGS += -DCPU_RECOMPUTE=0
else
  CXXFLAGS += -DCPU_RECOMPUTE=1
endif

ifeq ($(static),1)
  LIBRARY=libg2g.a
else
  LIBRARY=libg2g.so
endif

ifeq ($(full_double),1)
  CXXFLAGS += -DFULL_DOUBLE=1
  NVCCFLAGS += -DFULL_DOUBLE=1
else
  CXXFLAGS += -DFULL_DOUBLE=0
  NVCCFLAGS += -DFULL_DOUBLE=0
endif

## Define libraries
LIBRARIES := $(shell pkg-config gsl --libs) $(CUDA_LDFLAGS) $(CUDA_LIBS) -lrt

## Targets
all: $(LIBRARY)

include Makefile.depends

depend:
	@touch Makefile.depends
	@makedepend -Y -fMakefile.depends $(SRCS) > /dev/null 2>&1
	@sed -ri 's|(cuda/.*)\.o|\1.cu_o|g' Makefile.depends
	@rm -f Makefile.depends.bak

libg2g.a: $(OBJ)
	ar cru libg2g.a $(OBJ)
	ranlib libg2g.a

libg2g.so: $(OBJ)
# 	@echo "CUDAPATH:" $(CUDAPATH)
# 	@echo "CUDALINK:" $(CUDALINK)
	$(CXX) -shared $(LDFLAGS) -o libg2g.so $(OBJ) $(LIBRARIES)

clean:
# 	@if [ -f /opt/cuda/lib/libcudart.so ]; then \
# 	  echo "encontro"; else \
# 	  echo "no esta"; fi;
# 	@echo "CUDAPATH:" $(CUDAPATH)
# 	@echo "CUDALINK:" $(CUDALINK)
# 	make -C cuda clean
	rm -f *.o libg2g.so *.a cpu/*.o cuda/*.o
	rm -f cuda/*.cu_o cuda/*.cudafe* cuda/*.ptx cuda/*.hash cuda/*.cubin cuda/*.i cuda/*.ii cuda/*.fatbin.* cuda/*.cu.c
