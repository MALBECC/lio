################################################################################
# LIO MAKEFILE
################################################################################

all: liosolo liblio g2g


.PHONY: liosolo
liosolo: liblio
	$(MAKE) -C liosolo


.PHONY: liblio
liblio: g2g
	$(MAKE) -C lioamber


.PHONY: g2g
g2g:
	$(MAKE) -C g2g

.PHONY: hybrid
hybrid:
	$(MAKE) -C hybrid

.PHONY: clean
clean:
	$(MAKE) clean -C liosolo
	$(MAKE) clean -C lioamber
	$(MAKE) clean -C g2g
	$(MAKE) clean -C hybrid

################################################################################
