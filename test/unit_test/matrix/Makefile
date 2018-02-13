CFLAGS=-I. -I${HOME}/workspace/common/includes -I../../common/inc -I../../../g2g -I$(LIBXC_INCLUDES) -DGPU_KERNERLS=0
CLIBS=-L${HOME}/workspace/common/libs -lxc -lm -L../../../g2g -lg2g

CFLAGS_GPU=-I. -I${HOME}/workspace/common/includes -I../../common/inc -I../../../g2g -I$(LIBXC_INCLUDES) -DGPU_KERNELS=1 
CLIBS_GPU=-L${HOME}/workspace/common/libs -lxc -lm -L../../../g2g -lg2g

all: matrix_test

##matrix_test: matrix_test.cu
##	nvcc -G -g -o matrix_test matrix_test.cu $(CFLAGS_GPU) $(CLIBS) -DCPU_KERNELS

matrix_test: matrix_test.cu
	nvcc -G -g -o matrix_test matrix_test.cu $(CFLAGS_GPU) $(CLIBS_GPU) 

clean:
	rm -f *.o matrix_test

