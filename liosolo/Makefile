ifeq ($(gcc),1)
  CXX = gfortran
else
  CXX = ifort
endif 

 liosolo:
	$(CXX) main.f -L ../lioamber/ -llio-g2g -L ../g2g/ -lg2g -I../lioamber/ -o liosolo
 clean:
	rm liosolo

