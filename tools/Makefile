################################################################################

all: tdanalyze inputgen


.PHONY: inputgen
inputgen:
	$(MAKE) -C src/inputgen

.PHONY: tdanalyze
tdanalyze:
	$(MAKE) -C src/tdanalyze

.PHONY: clean
clean:
	$(MAKE) clean -C src/inputgen
	$(MAKE) clean -C src/tdanalyze

################################################################################
