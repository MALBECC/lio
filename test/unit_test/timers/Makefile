CFLAGS=-I. -I${HOME}/workspace/common/includes \
    -I../../common/inc \
    -I../../../g2g \
    -I$(LIBXC_INCLUDES) \
    -DGPU_KERNERLS=0
CLIBS=-L${HOME}/workspace/common/libs -lxc -lm -L../../../g2g -lg2g

CFLAGS_GPU=-I. -I${HOME}/workspace/common/includes \
    -I../../common/inc \
    -I../../../g2g \
    -I$(LIBXC_INCLUDES) \
    -DGPU_KERNELS=1

CLIBS_GPU=-L${HOME}/workspace/common/libs -lxc_cuda -lm -L../../../g2g -lg2g

all: timers_test

timers_test: timers_test.cu
	nvcc -G -g timers_test.cu $(CFLAGS_GPU) -o timers_test $(CLIBS_GPU)

clean:
	rm -f *.o timers_test

