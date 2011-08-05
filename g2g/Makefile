# -*- mode: make -*-
include Makefile.common
include cuda/Makefile.cuda

## Define sources and objects
SRCS:=$(wildcard *.cpp)
OBJ:=$(SRCS:%.cpp=%.o)

CUDA_LIBS:=-lcudart

ifeq ($(cpu),1)
  CPU_SOURCES  := $(wildcard cpu/*.cpp)
  SRCS         += $(CPU_SOURCES)
  OBJ          += $(CPU_SOURCES:%.cpp=%.o)
  CXXFLAGS     += -DCPU_KERNELS=1
else
  CUDA_SOURCES := $(wildcard cuda/*.cu)
  SRCS         += $(CUDA_SOURCES)
  OBJ          += $(CUDA_SOURCES:%.cu=%.cu_o)
  CUDA_LIBS         += -lcuda
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
	$(CXX) -shared $(LDFLAGS) -o libg2g.so $(OBJ) $(LIBRARIES)

clean:
	make -C cuda clean
	rm -f *.o *.so *.a cpu/*.o cuda/*.o
