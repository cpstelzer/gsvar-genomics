#!/bin/bash


OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

cd ..  # escape "SCRIPTS" folder

while read CLONE SAMPLENO FRAGSIZE

do # loop through clones/libraries

    echo "Currently processing  $CLONE"_"$SAMPLENO"	

    OUTMAIN="$CLONE/finalreads/"

    INGZ1="$OUTMAIN$CLONE"_"$SAMPLENO.roti-mito.R1.fq.gz"
    INGZ2="$OUTMAIN$CLONE"_"$SAMPLENO.roti-mito.R2.fq.gz"

    OUTREADS1="$OUTMAIN$CLONE"_"$SAMPLENO.roti-mito.R1.fq"
    OUTREADS2="$OUTMAIN$CLONE"_"$SAMPLENO.roti-mito.R2.fq"
    
    gunzip $INGZ1 $INGZ2
    
    ######### DO NOT CHANGE CODE BELOW! ####################################
    
    # Remove PCR-duplicates (Using FastUniq, i.e., mapped and unmapped fraction)
    module load fastuniq
    FU_READS1="$OUTMAIN$CLONE"_"$SAMPLENO.roti-mito.dedup.R1.fq"
    FU_READS2="$OUTMAIN$CLONE"_"$SAMPLENO.roti-mito.dedup.R2.fq"
    echo $OUTREADS1 > inlist
    echo $OUTREADS2 >> inlist
    fastuniq -i inlist -o $FU_READS1 -p $FU_READS2
    rm inlist

    # Quantify removal of reads/bases due to de-duplication
    FASTA_R1="temporary.R1.fasta"
    FASTA_R2="temporary.R2.fasta"
    COL[1]=$CLONE
    COL[2]=$SAMPLENO
    sed -n '1~4s/^@/>/p;2~4p' $OUTREADS1 > $FASTA_R1  # Convert to FASTA-format
    sed -n '1~4s/^@/>/p;2~4p' $OUTREADS2 > $FASTA_R2
    # Before duplicate removal 
    N_READS_R1=$(grep -c '^>' $FASTA_R1)  # number of reads in multi-fasta
    N_READS_R2=$(grep -c '^>' $FASTA_R1)
    (( COL[3] = $N_READS_R1 + $N_READS_R2 ))   # sum-up results of both mates
    TOTAL_BP_R1=$(grep -v '>' $FASTA_R1 | tr -d '\n' | wc -c) # total number of bases in fasta
    TOTAL_BP_R2=$(grep -v '>' $FASTA_R2 | tr -d '\n' | wc -c)
    (( COL[4] = $TOTAL_BP_R1 + $TOTAL_BP_R2 ))
    # After duplicate removal
    sed -n '1~4s/^@/>/p;2~4p' $FU_READS1 > $FASTA_R1  # Convert to FASTA-format
    sed -n '1~4s/^@/>/p;2~4p' $FU_READS2 > $FASTA_R2
    N_READS_R1=$(grep -c '^>' $FASTA_R1)  # number of reads in multi-fasta
    N_READS_R2=$(grep -c '^>' $FASTA_R1)
    (( COL[5] = $N_READS_R1 + $N_READS_R2 ))   # sum-up results of both mates
    TOTAL_BP_R1=$(grep -v '>' $FASTA_R1 | tr -d '\n' | wc -c) # total number of bases in fasta
    TOTAL_BP_R2=$(grep -v '>' $FASTA_R2 | tr -d '\n' | wc -c)
    (( COL[6] = $TOTAL_BP_R1 + $TOTAL_BP_R2 ))

    echo ${COL[*]} >> "SCRIPTS/fastuniq.test.csv"

    ############################# DO NOT CHANGE CODE ABOVE! ##############################

done < SCRIPTS/datasets/OHJ_rerun.csv


IFS=$OLDIFS

