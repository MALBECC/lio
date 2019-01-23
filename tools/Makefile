################################################################################

all: inputgen tdanalize


.PHONY: inputgen
inputgen:
	$(MAKE) -C inputgen

.PHONY: tdanalize
tdanalize:
	$(MAKE) -C tdanalize

.PHONY: clean
clean:
	$(MAKE) clean -C inputgen
	$(MAKE) clean -C tdanalize

################################################################################
