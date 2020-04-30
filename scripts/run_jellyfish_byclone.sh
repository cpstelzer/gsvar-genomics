#!/bin/bash

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

KMER=21                # kmer size
THREADS=8         # number of threads
MAX_RAM=1000000000        # maximum RAM
INMAIN="finalreads"
JELLYOUT="jellyfish_byclone"

module load jellyfish/2.1.4
cd ..

# Make a list with unique clone names
while read CLONE SAMPLENO FRAGSIZE
do
    echo "$CLONE" >> clones1.tmp
done < SCRIPTS/datasets/OHJ_all.csv
cat clones1.tmp | sort | uniq > clones2.tmp


while read CLONE
do
    
    echo "Currently processing  $CLONE"
    
    JELLYPATH="$CLONE"/"$INMAIN"/"$JELLYOUT"
    mkdir -p $JELLYPATH
 
    INFOLDER="$CLONE"/"$INMAIN"
    OUT_jf="$JELLYOUT"/"$CLONE"."$INMAIN".jf""
    OUT_hist="$JELLYOUT"/"$CLONE"."$INMAIN".histo""
        
    cd $INFOLDER

        echo "Prepare input FQ-files..."
        zcat *.fq.gz  > jelly_input.fq  # temporary input file
        
        echo "Run jellyfish kmer-counting..."
        jellyfish count -C -m $KMER -s $MAX_RAM -t $THREADS jelly_input.fq -o $OUT_jf
        jellyfish histo -t $THREADS $OUT_jf > $OUT_hist
        

        rm jelly_input.fq
    
    cd ../..


done < clones2.tmp

rm clones1.tmp clones2.tmp

IFS=$OLDIFS
