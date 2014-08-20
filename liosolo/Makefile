LIBS=-L/usr/lib -L/usr/lib64 -L../g2g -lg2g -lstdc++ \
     -L$(CUDA_HOME)/lib64/ -lcudart -lcublas \
     -L$(CUDA_HOME)/src/ -lfortran \

liosolo:
	ifort $(LIBS) -L ../lioamber -llio-g2g -L ../g2g -lg2g -I../lioamber main.f -o liosolo
clean:
	rm -rf liosolo
 

