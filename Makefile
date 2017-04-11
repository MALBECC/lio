################################################################################
# LIO MAKEFILE
################################################################################

all: liosolo/liosolo lioamber/liblio-g2g.so g2g/libg2g.so


.PHONY: liosolo
liosolo: liosolo/liosolo
liosolo/liosolo: lioamber/liblio-g2g.so
	$(MAKE) -C liosolo


.PHONY: liblio
liblio: lioamber/liblio-g2g.so
lioamber/liblio-g2g.so: g2g/libg2g.so
	$(MAKE) -C lioamber


.PHONY: g2g
g2g: g2g/libg2g.so
g2g/libg2g.so:
	$(MAKE) -C g2g


clean:
	$(MAKE) clean -C liosolo
	$(MAKE) clean -C lioamber
	$(MAKE) clean -C g2g

################################################################################
