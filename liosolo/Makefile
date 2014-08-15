FLAGS= -I../lioamber -L../g2g -L../lioamber
LIBRARIES=
ifeq ($(xeon_phi),1)
	FLAGS += -mmic
	LIBRARIES= -Qoption,link,-rpath=/opt/intel/composer_xe_2013_sp1.3.174/mkl/lib/mic -Qoption,link,-rpath=/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/mic/
endif

liosolo: main.f
	ifort $(FLAGS) -llio-g2g -lg2g main.f -o liosolo $(LIBRARIES)

clean:
	rm liosolo

