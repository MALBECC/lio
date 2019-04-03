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

#CFLAGS=-I. -I${HOME}/workspace/common/includes -I../../common/inc -I../../../g2g -I$(LIBXC_INCLUDES) -DGPU_KERNERLS=0
#CLIBS=-L${HOME}/workspace/common/libs -lxc -lm -L../../../g2g -lg2g

#CFLAGS_GPU=-I. -I${HOME}/workspace/common/includes -I../../common/inc -I../../../g2g -I$(LIBXC_INCLUDES) -DGPU_KERNELS=1 
#CLIBS_GPU=-L${HOME}/workspace/common/libs -lxc_cuda -lm -L../../../g2g -lg2g

ifeq ($(full_double),1)
    EXTRAS += -DFULL_DOUBLE=1
endif

all: main_gpu

main_cpu: main_cpu.cu
	nvcc -G -g -o main_cpu main_cpu.cu $(CFLAGS_GPU) $(CLIBS) -DCPU_KERNELS ${EXTRAS}

main_gpu: main_gpu.cu
	nvcc -G -g -o main_gpu main_gpu.cu $(CFLAGS_GPU) $(CLIBS_GPU) ${EXTRAS}

clean:
	rm -f *.o *.txt main_cpu main_gpu

