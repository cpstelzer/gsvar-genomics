#!/bin/bash

# Removing reads that map to the mitochondrial genome of B. plicatilis
# and keeps only reads that do NOT map. And subsequently, removing
# PCR duplicates (using FastUniq)

# C-P Stelzer 22.04.2019 (see page 74 of Notebook)

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

cd ..  # escape "SCRIPTS" folder
module load bedtools/2.27.1

MITOINDEX="mitoDNA/Bplic_mitogenome_IDX" # location of the mitochondrial Genome index

while read CLONE SAMPLENO FRAGSIZE

do # loop through clones/libraries

    echo "Currently processing  $CLONE"_"$SAMPLENO"	

    INMAIN="$CLONE/rotireads/"
    OUTMAIN="$CLONE/finalreads/"
    mkdir -p $OUTMAIN

    INREADS1_GZ="$INMAIN$CLONE"_"$SAMPLENO.rotireads.R1.fq.gz"
    INREADS2_GZ="$INMAIN$CLONE"_"$SAMPLENO.rotireads.R2.fq.gz"
    INREADS1="$INMAIN$CLONE"_"$SAMPLENO.rotireads.R1.fq"
    INREADS2="$INMAIN$CLONE"_"$SAMPLENO.rotireads.R2.fq"

    OUTREADS1="$OUTMAIN$CLONE"_"$SAMPLENO.roti-mito.R1.fq"
    OUTREADS2="$OUTMAIN$CLONE"_"$SAMPLENO.roti-mito.R2.fq"
    OUTLOG="$OUTMAIN$CLONE"_"$SAMPLENO.mitoDNA.bowtie2.log"        # bowtie2 logfile (contains info on %removed)

    echo "Unzipping read files..."
    gunzip $INREADS1_GZ $INREADS2_GZ
    
    # Declare names of intermediate files (not stored permanently!)
    OUTSAM="map2mito.tmp.sam"
    OUTBAM="map2mito.sorted.tmp.bam"
    OUTBAMI="$OUTBAM.bai"
    
    # Mapping to mitochondrial genome and filtering out reads that did NOT map
    echo "run bowtie2 ..."
    (bowtie2 -p 10 -X $FRAGSIZE -x $MITOINDEX -1 $INREADS1 -2 $INREADS2 -S $OUTSAM)2> $OUTLOG
    samtools view -S -h -@ 10 $OUTSAM | samtools sort -o $OUTBAM
    samtools index $OUTBAM

    samtools view -u -f 4 -F 264 $OUTBAM > tmps1.bam   # An unmapped read whose mate is mapped
    samtools view -u -f 8 -F 260 $OUTBAM > tmps2.bam   # A mapped read who’s mate is unmapped
    samtools view -u -f 12 -F 256 $OUTBAM > tmps3.bam  # Both reads of the pair are unmapped
    samtools merge unmapped.bam tmps1.bam tmps2.bam tmps3.bam
    samtools sort -n unmapped.bam -o unmapped.sorted.bam

    bedtools bamtofastq -i unmapped.sorted.bam -fq $OUTREADS1 -fq2 $OUTREADS2
    
    # Run jellyfish while reads are still uncompressed
    echo "Run jellyfish..."
	COMBIFQ="roti-mito.combined.tmp.fq"
	cat $OUTREADS1 $OUTREADS2 > $COMBIFQ
	JELLYPATH="$OUTMAIN"jellyfish/""
	mkdir -p $JELLYPATH	
	JELLYFILE="$JELLYPATH$CLONE"_"$SAMPLENO.jf"
	JELLYHISTO="$JELLYPATH$CLONE"_"$SAMPLENO.histo"
	jellyfish count -C -m 21 -s 1000000000 -t 10 $COMBIFQ -o $JELLYFILE
	jellyfish histo -t 1 $JELLYFILE > $JELLYHISTO

    # Run fastqc while reads are still uncompressed
    echo "Run fastqc..."
    FQCPATH="$OUTMAIN"fastqc/""
	mkdir -p $FQCPATH
	fastqc --threads 8 $OUTREADS1 $OUTREADS2 -o $FQCPATH >/dev/null

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
    
    echo ${COL[*]} >> "SCRIPTS/deduplication.report.csv"
    
    # Compression of read-files at the end of the run
    echo "Compressing read files..."    
    gzip $INREADS1 $INREADS2 
    gzip $OUTREADS1 $OUTREADS2
    gzip $FU_READS1 $FU_READS2

    # Remove all temporary files
    rm $OUTSAM $OUTBAM $OUTBAMI
    rm tmps1.bam tmps2.bam tmps3.bam
    rm unmapped.bam unmapped.sorted.bam
    rm $COMBIFQ
    
    echo "Completed for $CLONE"_"$SAMPLENO !"	


done < SCRIPTS/datasets/OHJ_all.csv


IFS=$OLDIFS


