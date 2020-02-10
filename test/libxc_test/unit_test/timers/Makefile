# Test for Libxc in the library folder
LIBXC = ${LIBXC_LIBS}
ifeq ("$(LIBXC)","")
    $(info >>> WARNING - Libxc not found in library folder. Please install Libxc before run the tests. <<<)
    HAS_LIBXC_INSTALLED := false
else
    HAS_LIBXC_INSTALLED := true
endif


CFLAGS=-I. -I${HOME}/workspace/common/includes \
    -I../../common/inc \
    -I../../../g2g \
    -I$(LIBXC_INCLUDES) \
    -DGPU_KERNERLS=0
CLIBS=-L${HOME}/workspace/common/libs -lxc -lm -L../../../g2g -lg2g

CFLAGS_GPU=-I. -I../../../g2g -I${LIBXC_INCLUDES} -DGPU_KERNELS=1
CLIBS_GPU=-L${LIBXC_LIBS} -lxc_cuda -lm -L../../../g2g -lg2g

all: timers_test

timers_test: timers_test.cu
	nvcc -G -g timers_test.cu $(CFLAGS_GPU) -o timers_test $(CLIBS_GPU)

clean:
	rm -f *.o timers_test

