################################################################################
# LIO MAKEFILE
################################################################################

all: liosolo liblio g2g tools


.PHONY: liosolo
liosolo: liblio
	$(MAKE) -C liosolo


.PHONY: liblio
liblio: g2g
	$(MAKE) -C lioamber


.PHONY: g2g
g2g:
	$(MAKE) -C g2g

.PHONY: tools
tools:
	$(MAKE) -C tools

.PHONY: check
check:
	$(MAKE) check -C test/tests_engine

.PHONY: clean
clean:
	$(MAKE) clean -C liosolo
	$(MAKE) clean -C lioamber
	$(MAKE) clean -C g2g
	$(MAKE) clean -C test/tests_engine
	$(MAKE) clean -C tools

################################################################################
