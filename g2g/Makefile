# -*- mode: make -*-
SRCS=$(wildcard *.cpp) $(wildcard cpu/*.cpp)
OBJ=$(SRCS:%.cpp=%.o)

CUDA_SOURCES=$(wildcard cuda/*.cu)
CUDA_HEADERS=$(wildcard cuda/*.h)

include cuda/Makefile.cuda
include Makefile.common

ifeq ($(cpu),1)
  CUDA_OBJ=
  CUDA_LIBS=
	CXXFLAGS += -DCPU_KERNELS=1
endif

ifeq ($(time), 1)
	CXXFLAGS += -DTIMINGS
endif

ifeq ($(histo), 1)
  CXXFLAGS += -DHISTOGRAM
endif

ifeq ($(profile),1)
	CXXFLAGS += -pg
	LDFLAGS += -pg
endif

ifeq ($(static),1)
	LIBRARY=libg2g.a
else
	LIBRARY=libg2g.so
endif

all: $(LIBRARY)

libg2g.a: $(CUDA_OBJ) $(OBJ)
	ar cru libg2g.a $(CUDA_OBJ) $(OBJ)
	ranlib libg2g.a
	
libg2g.so: $(CUDA_OBJ) $(OBJ)
	$(LINK) -shared $(LDFLAGS) $(EXTRA_FLAGS) -o libg2g.so $(OBJ) $(CUDA_OBJ) #potg.o

clean:
	make -C cuda clean
	rm -f *_test *.o *.so *.a *.cu_o cpu/*.o
