# Makefile for Harmonic_generation_py

CYAN := \033[1;36m
YELLOW := \033[1;33m
RED := \033[1;31m
RESET := \033[0m

PYTHONPATH := PYTHONPATH=$(shell pwd)/..
PARAMS := $(shell $(PYTHONPATH) python3 read_parameters.py)

PARAM_NAMES := EVOLVING_ATOM SAE_MODEL N_QN L_QN M_QN TOTAL_STATES N R_MAX L_MAP CONFINED CONF_MODEL L_MAX K_MAX DT TIME_STEP STATE_SYMB L
$(foreach i,1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17,$(eval $(word $i,$(PARAM_NAMES)) := $(word $i,$(PARAMS))))


COLLOC_dPATH := collocation_points/generator
INPUT_dPATH  := GPSM_states_S-matrix/data_GPSM_states_S-matrix
TEVO_DAT     := Time_evolution_data

COLLOC_DAT := $(COLLOC_dPATH)/Algo-3_N=$(N)_Gauss_Lobatto_collocation_points.txt

ifeq (${CONFINED}, False)
    GPSM_DAT := $(INPUT_dPATH)/Free_atom/${EVOLVING_ATOM}_States_${SAE_MODEL}_l=${L_QN}_nos=${TOTAL_STATES}_N=${N}_rmax=${R_MAX}_Lmap=${L_MAP}.dat
    SMAT_DAT := $(INPUT_dPATH)/Free_atom/${EVOLVING_ATOM}_Smatrix_${SAE_MODEL}_m=${M_QN}_lmax=${L_MAX}_kmax=${K_MAX}_N=${N}_rmax=${R_MAX}_Lmap=${L_MAP}_dt=${DT}.npy
    TEVO_DAT := $(TEVO_DAT)/VEvo_nopt=${TIME_STEP}_${EVOLVING_ATOM}(${STATE_SYMB})_m=${M_QN}_${SAE_MODEL}_L=${L}_kmax=${K_MAX}_N=${N}_rmax=${R_MAX}_Lmap=${L_MAP}_dt=${DT}.dat
else
    GPSM_DAT := $(INPUT_dPATH)/Confined_atom/${EVOLVING_ATOM}@C60_States_${SAE_MODEL}_${CONF_MODEL}_l=${L_QN}_nos=${TOTAL_STATES}_N=${N}_rmax=${R_MAX}_Lmap=${L_MAP}.dat
    SMAT_DAT := $(INPUT_dPATH)/Confined_atom/${EVOLVING_ATOM}@C60_Smatrix_${SAE_MODEL}_${CONF_MODEL}_m=${M_QN}_lmax=${L_MAX}_kmax=${K_MAX}_N=${N}_rmax=${R_MAX}_Lmap=${L_MAP}_dt=${DT}.npy
    TEVO_DAT := $(TEVO_DAT)/VEvo_nopt=${TIME_STEP}_${EVOLVING_ATOM}(${STATE_SYMB})@C60_m=${M_QN}_${SAE_MODEL}_${CONF_MODEL}_L=${L}_kmax=${K_MAX}_N=${N}_rmax=${R_MAX}_Lmap=${L_MAP}_dt=${DT}.dat
endif


.PHONY : all HHG_spectra clear_colloc_pt clear_GPSM_SMAT clear_evo_data clobber

all : $(COLLOC_DAT) $(GPSM_DAT) $(SMAT_DAT) $(TEVO_DAT) HHG_spectra


$(COLLOC_DAT) :
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|               Algo-3_Gauss_Lobatto_colloc_pt.py                  |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)python3 $(COLLOC_dPATH)/Algo-3_Gauss_Lobatto_colloc_pt.py -N $(N)$(RESET)"
	@python3 $(COLLOC_dPATH)/Algo-3_Gauss_Lobatto_colloc_pt.py -N $(N)

$(GPSM_DAT) :
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                    GPSM_state_generator.py                       |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)python3 GPSM_states_S-matrix/GPSM_states_generator.py -v$(RESET)"
	@python3 GPSM_states_S-matrix/GPSM_states_generator.py -v

$(SMAT_DAT) :
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                     S_matrix_generator.py                        |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)python3 GPSM_states_S-matrix/S_matrix_generator.py$(RESET)"
	@python3 GPSM_states_S-matrix/S_matrix_generator.py

$(TEVO_DAT) :
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                    vector_time_evolution.py                      |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)python3 vector_time_evolution.py -A$(RESET)"
	@python3 vector_time_evolution.py -A

HHG_spectra :
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(CYAN)|                          HHG_spectra.py                          |$(RESET)"
	@echo "$(CYAN)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$(RESET)"
	@echo "$(YELLOW)python3 HHG_spectra_analysis/HHG_spectra.py -A$(RESET)"
	@python3 HHG_spectra_analysis/HHG_spectra.py -A

clear_colloc_pt :
	@echo "$(RED)rm -f collocation_points/generator/*.txt$(RESET)"
	@rm -f collocation_points/generator/*.txt

clear_GPSM_SMAT :
	@echo "$(RED)rm -f GPSM_states_S-matrix/data_GPSM_states_S-matrix/Free_atom/*$(RESET)"
	@rm -f GPSM_states_S-matrix/data_GPSM_states_S-matrix/Free_atom/*
	@echo "$(RED)rm -f GPSM_states_S-matrix/data_GPSM_states_S-matrix/Confined_atom/*$(RESET)"
	@rm -f GPSM_states_S-matrix/data_GPSM_states_S-matrix/Confined_atom/*

clear_evo_data :
	@echo "$(RED)rm -f Time_evolution_data/*$(RESET)"
	@rm -f Time_evolution_data/*

clobber : clear_colloc_pt clear_GPSM_SMAT clear_evo_data
	@echo "$(RED)All generated data cleared!$(RESET)"
