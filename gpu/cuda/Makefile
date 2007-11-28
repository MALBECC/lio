# -*- mode: make -*-

CUDA_SOURCES=$(wildcard *.cu)
include Makefile.cuda

all: $(CUDA_OBJ)

clean:
	rm -f *.cu_o *.cudafe* *.ptx *.hash *.cubin *.i *.ii *.fatbin.* *.cu.c
