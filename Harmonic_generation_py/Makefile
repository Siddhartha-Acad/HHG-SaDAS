# Makefile for Harmonic_generation_py

PYTHONPATH   := PYTHONPATH=$(shell pwd)/..
COLLOC_dPATH := collocation_points/generator
GPSM_dPATH   := GPSM_states_S-matrix/data_GPSM_states_S-matrix/Free_atom

PARAMS := $(shell $(PYTHONPATH) python3 read_parameters.py)

PARAM_NAMES := EVOLVING_ATOM SAE_MODEL L TOTAL_STATES N R_MAX L_MAP CONFINED
${foreach i,1 2 3 4 5 6 7 8,${eval ${word $i,${PARAM_NAMES}} := ${word $i,${PARAMS}}}}


COLLOC_DAT := $(COLLOC_dPATH)/Algo-3_N=$(N)_Gauss_Lobatto_collocation_points.txt

ifeq (${CONFINED}, False)
    GPSM_DAT := $(GPSM_dPATH)/${EVOLVING_ATOM}_States_${SAE_MODEL}_l=${L}_nos=${TOTAL_STATES}_N=${N}_rmax=${R_MAX}_Lmap=${L_MAP}.dat
endif


.PHONY : all clear_colloc_pt

all : $(COLLOC_DAT) $(GPSM_DAT)

$(COLLOC_DAT) :
	python3 $(COLLOC_dPATH)/Algo-3_Gauss_Lobatto_colloc_pt.py -N $(N)

$(GPSM_DAT) :
	python3 GPSM_states_S-matrix/GPSM_states_generator.py

clear_colloc_pt :
	rm -f collocation_points/generator/*.txt

clear_GPSM :
	rm -f GPSM_states_S-matrix/data_GPSM_states_S-matrix/Free_atom/*States*.dat

