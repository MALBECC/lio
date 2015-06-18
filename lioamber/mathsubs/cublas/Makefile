include ../../make.inc

all: $(obj_path)fortran.o 

$(obj_path)fortran.o: $(CUDA_HOME)/src/fortran.c
	$(CC) -fPIC -DCUBLAS_GFORTRAN -O3 -c $< -o $@ -I$(CUDA_HOME)/include

#$(obj_path)cublasmath.o: $(obj_path)garcha_mod.mod $(obj_path)fortran.o
#	$(FC) -module ../../obj  -fpp -Dpack -fPIC -DG2G -g -c -O3 $(*F).f -o $@


clean:
	rm -f $(obj_path)fortran.o;
