# Test for Libxc in the library folder
LIBXC = ${LIBXC_LIBS}
ifeq ("$(LIBXC)","")
    $(info >>> WARNING - Libxc not found in library folder. Please install Libxc before run the tests. <<<)
    HAS_LIBXC_INSTALLED := false
else
    HAS_LIBXC_INSTALLED := true
endif

CFLAGS=-I${LIBXC_INCLUDES} -I../commons -DGPU_KERNERLS=0
CLIBS=-L${LIBXC_LIBS} -lxc -lm -L../../../g2g -lg2g

CFLAGS_GPU=-I${LIBXC_INCLUDES} -I../commons -I../../../g2g -DGPU_KERNELS=1 --ptxas-options=-v
CLIBS_GPU=-L${LIBXC_LIBS} -lxc_cuda -lm -L../../../g2g -lg2g

all: matrix_test

matrix_test: matrix_test.cu
	nvcc -G -g matrix_test.cu $(CFLAGS_GPU) -o matrix_test $(CLIBS_GPU)

clean:
	rm -f *.o matrix_test

