# Makefile for Harmonic_generation_py

PROJECT_ROOT := $(shell pwd)/..
COLLOC_SCRIPT := collocation_points/generator

N := $(shell PYTHONPATH=$(PROJECT_ROOT) python3 -c "from parameters import N; print(N)")

COLLOC_DAT := $(COLLOC_SCRIPT)/Algo-3_N=$(N)_Gauss_Lobatto_collocation_points.txt

.PHONY: clear_colloc_pt


$(COLLOC_DAT):
	python3 $(COLLOC_SCRIPT)/Algo-3_Gauss_Lobatto_colloc_pt.py -N $(N)

clear_colloc_pt:
	rm -f collocation_points/generator/*.txt

