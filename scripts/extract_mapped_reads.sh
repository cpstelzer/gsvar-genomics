#!/bin/bash

# Extract mapped reads from BAM alignment file and store them in folder $CLONE/aqtrim/mapped

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

cd ..  # escape "SCRIPTS" folder
module load bedtools/2.27.1

while read CLONE SAMPLENO FRAGSIZE

do

	echo "Currently processing  $CLONE"	

    MAPPATH="$CLONE/aqtrim/mapped/"
    
    mkdir -p $MAPPATH

    BAMFILE="$CLONE"/mapping/"$CLONE"_"$SAMPLENO"_to_VBCFpol.sorted.bam"" # the original alignment to the reference
    samtools view -u -F4 $BAMFILE > TEMP_mapped.bam # extract only reads that map to the reference genome 
    samtools sort -n TEMP_mapped.bam -o TEMP_mapped.sorted.bam
    FQ_R1="$MAPPATH$CLONE"_"$SAMPLENO".mapped.R1.fq"" # uncompressed fq
    FQ_R2="$MAPPATH$CLONE"_"$SAMPLENO".mapped.R2.fq""
    bedtools bamtofastq -i TEMP_mapped.sorted.bam -fq $FQ_R1  -fq2 $FQ_R2 2>/dev/null
    gzip $FQ_R1 $FQ_R2

    rm TEMP_mapped.bam TEMP_mapped.sorted.bam
    
done < SCRIPTS/datasets/OHJ7_mbl2019.csv


IFS=$OLDIFS
