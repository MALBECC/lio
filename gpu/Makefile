# -*- mode: make -*-

OBJ=init.o timer.o matrix.o partition.o functions.o weight.o

CUDA_SOURCES=$(wildcard cuda/*.cu)
CUDA_HEADERS=$(wildcard cuda/*.h)

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

all: $(LIBRARY) #accum_test dist_test exchnum_test

#%_test: %_test.o libg2g.so 
#	$(LINK) $(LDFLAGS) -o $@ $(OBJ) $(LIB) $@.o -lg2g -L.

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

