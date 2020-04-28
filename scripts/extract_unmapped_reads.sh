#!/bin/bash

# Code based on: https://github.com/alvaralmstedt/Tutorials/wiki/Separating-mapped-and-unmapped-reads-from-libraries
# See notebook page 37 (10.3.2020)

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

cd ..  # escape SCRIPTS folder

ASSEMBLY="VBCFpol"   # Polished genome. Alternatively, VBCF (=unpolished version)

module load bedtools/2.27.1

while read CLONE SAMPLENO FRAGSIZE

do
    echo "Currently processing  $CLONE"	
	
    MAPPATH="$CLONE/mapping/"
    TRIMMPATH="$CLONE/aqtrim/"

    mkdir -p "$TRIMMPATH/unmapped/"

    INFILE="$MAPPATH$CLONE"_"$SAMPLENO"_to_"$ASSEMBLY.sorted.bam"
	
    # Filtering out reads that did not map (code from novocraft.com)
    samtools view -u -f 4 -F 264 $INFILE > tmps1.bam       # An unmapped read whose mate is mapped
    samtools view -u -f 8 -F 260 $INFILE > tmps2.bam       # A mapped read who’s mate is unmapped
    samtools view -u -f 12 -F 256 $INFILE > tmps3.bam      # Both reads of the pair are unmapped

    samtools merge unmapped.bam tmps1.bam tmps2.bam tmps3.bam
    samtools sort -n unmapped.bam -o unmapped.sorted.bam

    bedtools bamtofastq -i unmapped.sorted.bam -fq unmapped.R1.fq  -fq2 unmapped.R2.fq

    TRIMMED_R1="$TRIMMPATH$CLONE"_"$SAMPLENO.trimmed.R1.fq.gz"
    TRIMMED_R2="$TRIMMPATH$CLONE"_"$SAMPLENO.trimmed.R2.fq.gz"
    UNMAPPED_R1="$TRIMMPATH"/unmapped/"$CLONE"_"$SAMPLENO.unmapped.R1.fq"
    UNMAPPED_R2="$TRIMMPATH"/unmapped/"$CLONE"_"$SAMPLENO.unmapped.R2.fq"
	
    bedtools bamtofastq -i unmapped.sorted.bam -fq $UNMAPPED_R1  -fq2 $UNMAPPED_R2
 
    # Note: all nonmapped reads should be also stored in a central folder!!
	
    # Creating a report of the unmapped reads
    COL[1]=$CLONE
    COL[2]=$SAMPLENO 
    COL[3]=$(zcat $TRIMMED_R1 | grep -c "+")
    COL[4]=$(zcat $TRIMMED_R2 | grep -c "+")
    COL[5]=$(grep -c "+" $UNMAPPED_R1)
    COL[6]=$(grep -c "+" $UNMAPPED_R2)
    
    echo ${COL[*]} 
    IFS=$OLDIFS
    echo ${COL[*]} >> "SCRIPTS/unmapped.reads.csv"
    IFS=$NEWIFS

    # Delete and zip intermediate files to save space
    rm tmps*
    rm unmapped.bam
    rm unmapped.sorted.bam
    gzip $UNMAPPED_R1
    gzip $UNMAPPED_R2


done < SCRIPTS/datasets/OHJ7_mbl2019.csv


IFS=$OLDIFS
