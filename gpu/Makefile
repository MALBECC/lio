# -*- mode: make -*-

OBJ=timer.o
CUDA_OBJ=$(wildcard cuda/*.cu_o)

include cuda/Makefile.cuda
include Makefile.common

all: libg2g.so

cuda_test: libg2g.so main.o
	$(LINK) $(LDFLAGS) -o cuda_test $(OBJ) $(LIB) main.o -lg2g -L.

libg2g.so: $(OBJ) $(CUDA_OBJ)
	$(MAKE) -C cuda
	$(LINK) -shared $(LDFLAGS) -o libg2g.so $(OBJ) $(CUDA_OBJ) $(LIB)

clean:
	make -C cuda clean
	rm -f *.o *.so

