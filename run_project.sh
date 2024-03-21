#!/bin/bash
# script for running variant clustering project

# Flags for which files to run
CALCULATE_POSITIONS=1	# whether to run step 1
RUN_SIMULATIONS=1	# whether to run simulations
EXTRACT_CLINVAR=1	# whether to process clinvar variants to compare to UKBB ones
RUN_CLUSTERSING=1	# whether to run clustering step

# load in file locations
. file_locations
# create output directories
mkdir -p results
mkdir -p plots

loadR

# step 1
# calculate relative positions within gene/transcript for each variant
if [[ $CALCULATE_POSITIONS ]]
then
	for i in {1..22} X	# first collate the annotations for all variants
	do 
		python code/Step1_calculate_positions/parse_variant_consequences.py --chr $i
	done
	Rscript ./code/Step1_calculate_positions/calculate_variant_positions.R # then calculate relative positions of variants
fi
if [[ $RUN_SIMULATIONS ]]
then
	Rscript code/Step2_simulations/simulate_genes.R
fi
if [[ $EXTRACT_CLINVAR ]]
then
	# get pathogenic variants from Clinvar VCF files
	python code/Step3_extract_clinvar/extract_path_variants.py
	# process these to get relative positions and then cluster
	Rscript code/Step3_extract_clinvar/calculate_clinvar_relative_positions.R
	python code/Step3_extract_clinvar/extract_possible_variants.R
fi
if [[ $RUN_CLUSTERING ]]
then
	python ./code/Step4_cluster_genes/cluster_genes.py
fi