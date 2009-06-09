# -*- mode: make -*-

OBJ=init.o timer.o matrix.o partition.o functions.o weight.o

CUDA_SOURCES=$(wildcard cuda/*.cu)
CUDA_HEADERS=$(wildcard cuda/*.h)

ifeq ($(cpu),1)
	CXXFLAGS += -DCPU_KERNELS=1
endif

include cuda/Makefile.cuda
include Makefile.common

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

ifeq ($(store),0)
	CXXFLAGS += -DSTORE_FUNCTIONS=0
else
	ifeq ($(store),1)
		CXXFLAGS += -DSTORE_FUNCTIONS=1 -DDENSITY_ALT=0
	else
		CXXFLAGS += -DSTORE_FUNCTIONS=1 -DDENSITY_ALT=1
	endif
endif

all: $(LIBRARY)

libg2g.a: $(CUDA_OBJ) $(OBJ)
	$(MAKE) -C cuda
	ar cru libg2g.a $(CUDA_OBJ) $(OBJ)
	ranlib libg2g.a
	
libg2g.so: $(CUDA_OBJ) $(OBJ)
	$(MAKE) -C cuda
	$(LINK) -shared $(LDFLAGS) $(EXTRA_FLAGS) -o libg2g.so $(OBJ) $(CUDA_OBJ) $(LIB)

clean:
	make -C cuda clean
	rm -f *_test *.o *.o *.so *.a *.cu_o

