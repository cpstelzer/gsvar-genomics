#!/bin/bash

# Verwenden Sie diese Datei so (via submit_hecaton_job.sge)
# cd $SCRATCH/hecaton_singularity
# ./hecaton.sh bash hecaton_preprocess.sh

bash /hecaton/bash/preprocess.sh $SCRATCH/data/assembly/OHJ7i3n10_VBCFpol_sorted.fa
