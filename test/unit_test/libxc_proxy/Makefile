CC=gcc
CFLAGS=-I. -I${HOME}/workspace/common/includes -I../../common/inc -I../../../g2g -I$(LIBXC_INCLUDES)
CLIBS=-L${HOME}/workspace/common/libs -lxc -lm -L../../../g2g -lg2g

NVCC=nvcc
CFLAGS_GPU=-I. -I${HOME}/workspace/common/includes -I../../common/inc -I../../../g2g -I$(LIBXC_INCLUDES) -DGPU_KERNELS=1 --ptxas-options=-v
CLIBS_GPU=-L${HOME}/workspace/common/libs -lxc_cuda -lm -L../../../g2g -lg2g

ifeq ($(full_double),1)
    EXTRAS += -DFULL_DOUBLE=1
endif

all:  main_gpu

main_cpu: main_cpu.cpp 
	nvcc -G -g -o main_cpu main_cpu.cpp $(CFLAGS) $(CLIBS) -DCPU_KERNELS

main_gpu: main_gpu.cu
	nvcc -G -g -o main_gpu main_gpu.cu $(CFLAGS_GPU) $(CLIBS_GPU) $(EXTRAS)

clean:
	rm -f *.o main_cpu main_gpu

