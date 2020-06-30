#!/bin/bash

OLDIFS=$IFS     # save the existing field separator
IFS=";"

cd ..  # escape SCRIPTS folder


INMAIN="/finalreads2/"
OUTMAIN="/finalreads3/"

module load bbtools/38.34

while read CLONE SAMPLENO FRAGSIZE

do
	echo "Currently processing $CLONE $SAMPLENO"	

    INPATH="$CLONE$INMAIN"
    OUTPATH="$CLONE$OUTMAIN"
    mkdir -p $OUTPATH

    INGZ1="$INPATH$CLONE"_"$SAMPLENO".finalreads.R1.fq.gz""
    INGZ2="$INPATH$CLONE"_"$SAMPLENO".finalreads.R2.fq.gz""

    OUTGZ1="$OUTPATH$CLONE"_"$SAMPLENO".finalreads.R1.fq.gz""
    OUTGZ2="$OUTPATH$CLONE"_"$SAMPLENO".finalreads.R2.fq.gz""

    OUTREPORT="$OUTPATH$CLONE"_"$SAMPLENO".filtertiles.report.txt""

    # filterbytile.sh in1=$INGZ1 in2=$INGZ2 out1=$OUTGZ1 out2=$OUTGZ2 &> $OUTREPORT
    # Using more aggressive filtering option (https://www.biostars.org/p/228762/)
    filterbytile.sh in1=$INGZ1 in2=$INGZ2 out1=$OUTGZ1 out2=$OUTGZ2 ud=0.75 qd=1 ed=1 ua=.5 qa=.5 ea=.5 &> $OUTREPORT

bbduk.sh in=reads.fq out=clean.fq qtrim=r trimq=10


    # Run FASTQC and JELLYFISH on filtered reads

    OUTFQ1="$OUTPATH$CLONE"_"$SAMPLENO".finalreads.R1.fq""
    OUTFQ2="$OUTPATH$CLONE"_"$SAMPLENO".finalreads.R2.fq""
    gunzip $OUTGZ1 $OUTGZ2

    FQCPATH="$OUTPATH"fastqc/""
    mkdir -p $FQCPATH
    fastqc --threads 8 $OUTFQ1 $OUTFQ2 -o $FQCPATH

    JELLYPATH="$OUTPATH"jellyfish""
	mkdir -p $JELLYPATH	
	JELLYFILE="$JELLYPATH/$CLONE"_"$SAMPLENO.jf"
	JELLYHISTO="$JELLYPATH/$CLONE"_"$SAMPLENO.histo"
	jellyfish count -C -m 21 -s 1000000000 -t 8 $COMBIFQ -o $JELLYFILE
	jellyfish histo -t 1 $JELLYFILE > $JELLYHISTO

    # Count reads and BPs in filtered reads

    FASTA_R1="temp.R1.fasta"
    FASTA_R2="temp.R2.fasta"
    COL[1]=$CLONE
    COL[2]=$SAMPLENO
    sed -n '1~4s/^@/>/p;2~4p' $OUTFQ1 > $FASTA_R1  # Convert to FASTA-format
    sed -n '1~4s/^@/>/p;2~4p' $OUTFQ2 > $FASTA_R2
    N_READS_R1=$(grep -c '^>' $FASTA_R1)  # number of reads in multi-fasta
    N_READS_R2=$(grep -c '^>' $FASTA_R2)
    (( COL[3] = $N_READS_R1 + $N_READS_R2 ))   # sum-up results of both mates
    TOTAL_BP_R1=$(grep -v '>' $FASTA_R1 | tr -d '\n' | wc -c) # total number of bases in fasta
    TOTAL_BP_R2=$(grep -v '>' $FASTA_R2 | tr -d '\n' | wc -c)
    (( COL[4] = $TOTAL_BP_R1 + $TOTAL_BP_R2 ))
    echo ${COL[*]} >> "SCRIPTS/finalreads3.report.csv"


    gzip $OUTFQ1 $OUTFQ2
    rm temp*

done < SCRIPTS/datasets/OHJ_wbadtiles.csv



IFS=$OLDIFS
