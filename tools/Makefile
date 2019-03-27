################################################################################

all: tdanalyze inputgen bdens

.PHONY: inputgen
inputgen:
	$(MAKE) -C src/inputgen

.PHONY: tdanalyze
tdanalyze:
	$(MAKE) -C src/tdanalyze

.PHONY: bdens
bdens:
	$(MAKE) -C src/buildens

.PHONY: clean
clean:
	$(MAKE) clean -C src/buildens
	$(MAKE) clean -C src/inputgen
	$(MAKE) clean -C src/tdanalyze

################################################################################
