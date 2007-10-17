# -*- mode: make -*-

OBJ=timer.o matrix.o common.o

CUDA_SOURCES=$(wildcard cuda/*.cu)

include cuda/Makefile.cuda
include Makefile.common

all: libg2g.so dist_test accum_test

%_test: libg2g.so %_test.o
	$(LINK) $(LDFLAGS) -o $@ $(OBJ) $(LIB) $@.o -lg2g -L.

libg2g.so: $(OBJ) $(CUDA_OBJ)
	$(MAKE) -C cuda
	$(LINK) -shared $(LDFLAGS) -o libg2g.so $(OBJ) $(CUDA_OBJ) $(LIB)

clean:
	make -C cuda clean
	rm -f tests/*_test tests/*.o *.o *.so

