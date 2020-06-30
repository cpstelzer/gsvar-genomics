#!/bin/bash

LIBNAME="IK1_vbcf73392"
OUTDIR="$SCRATCH/data/hecaton_output/$LIBNAME"
mkdir -p $OUTDIR
cd /hecaton
nextflow run -c nextflow/nextflow.config -w $SCRATCH/hecaton_singularity nextflow/hecaton.nf --genome_file $SCRATCH/data/assembly/OHJ7i3n10_VBCFpol_sorted.fa --reads "$SCRATCH/data/finalreads/$LIBNAME.roti-mito.dedup.R{1,2}.fq.gz" --manta_config docker/configManta_weight_1.py.ini --output_dir $OUTDIR --model_file models/random_forest_model_concat_A_thaliana_ColxCvi_O_sativa_Suijing18_coverage_10x_insertions_balanced_subsample.pk
