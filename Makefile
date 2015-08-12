.PHONY: all clean hooks

ifeq ($(cublas),1)
all:
#	cd g2g; make;
	cd lioamber/; mkdir -p obj;
	cd lioamber/mathsubs/cublas; make;
	cd lioamber; make;
	cd liosolo; make;

clean:
#	cd g2g; make clean;
	cd lioamber/mathsubs/cublas; make clean;
	cd lioamber; make clean;
	cd liosolo; make clean;
else
all:
#	cd g2g; make;
	cd lioamber; make;
	cd liosolo; make;

clean:
#	cd g2g; make clean;
	cd lioamber; make clean;
	cd liosolo; make clean;
endif
