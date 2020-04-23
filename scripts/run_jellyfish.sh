#!/bin/bash

# Run jellyfish with kmer-size of 21 on "rotifer reads"
# See notebook, page 65 (8.4.2020)

# "Rotifer reads" are all mapped reads, plus all unmapped reads that
# have not been identified as contaminants

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

KMER=21                # kmer size
THREADS=6         # number of threads
MAX_RAM=1000000000        # maximum RAM

INMAIN="/rotireads/"

cd ..  # escape "SCRIPTS/" folder
module load jellyfish/2.1.4

while read CLONE SAMPLENO FRAGSIZE

do

        echo "Currently processing  $CLONE"

    JELLYPATH="$CLONE$INMAIN"jellyfish/""
    mkdir -p $JELLYPATH

    IN1="$CLONE$INMAIN$CLONE"_"$SAMPLENO".rotireads.R1.fq.gz""
    IN2="$CLONE$INMAIN$CLONE"_"$SAMPLENO".rotireads.R2.fq.gz""
    OUT_jf="$JELLYPATH$CLONE"_"$SAMPLENO".rotireads.jf""
    OUT_hist="$JELLYPATH$CLONE"_"$SAMPLENO".rotireads.histo""

    echo "Prepare input FQ-file..."
    zcat $IN1 $IN2  > jelly_input.fq  # temporary input file

    echo "Run jellyfish kmer-counting..."
    jellyfish count -C -m $KMER -s $MAX_RAM -t $THREADS jelly_input.fq -o $OUT_jf
    jellyfish histo -t $THREADS $OUT_jf > $OUT_hist

    rm jelly_input.fq


done < SCRIPTS/datasets/OHJ_all.csv


IFS=$OLDIFS


