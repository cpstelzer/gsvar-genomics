#!/bin/bash

# Run jellyfish with kmer-size of 21 on "rotifer reads"
# and fastqc on a small selection of clones/libraries
# See notebook, page 76 (23.4.2020)


OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

KMER=21                # kmer size
THREADS=8         # number of threads
MAX_RAM=1000000000        # maximum RAM

INMAIN="/finalreads/"

cd ..  # escape "SCRIPTS/" folder
module load jellyfish/2.1.4

while read CLONE SAMPLENO FRAGSIZE

do

        echo "Currently processing  $CLONE"

    JELLYPATH="$CLONE$INMAIN"jellyfish/""
    mkdir -p $JELLYPATH

    IN1_GZ="$CLONE$INMAIN$CLONE"_"$SAMPLENO".roti-mito.dedup.R1.fq.gz""
    IN2_GZ="$CLONE$INMAIN$CLONE"_"$SAMPLENO".roti-mito.dedup.R2.fq.gz""
    IN1="$CLONE$INMAIN$CLONE"_"$SAMPLENO".roti-mito.dedup.R1.fq""
    IN2="$CLONE$INMAIN$CLONE"_"$SAMPLENO".roti-mito.dedup.R2.fq""
    OUT_jf="$JELLYPATH$CLONE"_"$SAMPLENO".roti-mito.dedup.jf""
    OUT_hist="$JELLYPATH$CLONE"_"$SAMPLENO".rotir-mito.dedup.histo""

    echo "Prepare input FQ-file..."
    gunzip $IN1_GZ $IN2_GZ

    cat $IN1 $IN2  > jelly_input.fq  # temporary input file

    echo "Run jellyfish kmer-counting..."
    jellyfish count -C -m $KMER -s $MAX_RAM -t $THREADS jelly_input.fq -o $OUT_jf
    jellyfish histo -t $THREADS $OUT_jf > $OUT_hist

    rm jelly_input.fq

    echo "Run fastqc..."
    FQCPATH="$INMAIN"fastqc/""
	mkdir -p $FQCPATH
	fastqc --threads 8 $IN1 $IN2 -o $FQCPATH >/dev/null

    gzip $IN1 $IN2


done < SCRIPTS/datasets/OHJ_IK1_OHJ104.csv


IFS=$OLDIFS


