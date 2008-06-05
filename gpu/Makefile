# -*- mode: make -*-

OBJ=timer.o matrix.o common.o

CUDA_SOURCES=$(wildcard cuda/*.cu)
CUDA_HEADERS=$(wildcard cuda/*.h)

include cuda/Makefile.cuda
include Makefile.common

ifeq ($(time), 1)
	CXXFLAGS += -DDO_TIMINGS
endif


all: libg2g.so #accum_test dist_test exchnum_test

%_test: %_test.o libg2g.so
	$(LINK) $(LDFLAGS) -o $@ $(OBJ) $(LIB) $@.o -lg2g -L.

# $(CUDA_HEADERS)
libg2g.so: $(CUDA_OBJ) $(OBJ) 
	$(MAKE) -C cuda
	$(LINK) -shared $(LDFLAGS) $(EXTRA_FLAGS) -o libg2g.so $(OBJ) $(CUDA_OBJ) $(LIB)

clean:
	make -C cuda clean
	rm -f *_test *.o *.o *.so *.cu_o *.i *.ii *.cudafe?.* *.cubin *.fatbin.* *.hash *.ptx *.cu.cpp

