#!/bin/bash

# Subtract Tetraselmis-contamination (i.e., unmapped reads of finalmapping that DID MAP to the Tetraselmis striata genome) from finalreads
# 1) Load Unmapped_finalreads_to_Tetra_alignment (BAM-file)
# 2) Extract Tetra-unmapped reads (these are, by definition rotifer reads)
# 3) Merge these "rotifer reads" with mapped reads of the final assembly
# 4) Run fastqc and jellyfish k-mer analysis on final reads
# 5) Store the joint set as 'finalreads2' (overwrite previous files, which include Tetra-contamination; for safety, I still have the previous versions
#    as zipped archives finalreads & finalmapping.

# Written before, but likely to be run on 27.05.2020


OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

cd ..  # escape contaminants folder

module load bedtools/2.27.1

while read CLONE SAMPLENO FRAGSIZE

do
    echo "Currently processing  $CLONE"_" $SAMPLENO"	
	
    OUTMAIN="$CLONE/finalreads2/"
    mkdir -p $OUTMAIN   
    
    # 1) Load Unmapped_finalreads_to_Tetra_alignment (BAM-file)
    TETRAPATH="$CLONE/finalmapping/tetraselmis/"
    INFILE="$TETRAPATH$CLONE"_"$SAMPLENO".unmapped_finalreads_2_tetra.sorted.bam""  # Unmapped finalreads aligned to T. striata

    # 2) Extract Tetra-unmapped reads (these are, by definition, rotifer reads)
    samtools view -u -f 4 -F 264 $INFILE > tffs1.bam       # An unmapped read whose mate is mapped
    samtools view -u -f 8 -F 260 $INFILE > tffs2.bam       # A mapped read who’s mate is unmapped
    samtools view -u -f 12 -F 256 $INFILE > tffs3.bam      # Both reads of the pair are unmapped
    samtools merge TR_TF.unmapped.bam tffs1.bam tffs2.bam tffs3.bam
    samtools sort -n TR_TF.unmapped.bam -o TR_TF.unmapped.sorted.bam
    bedtools bamtofastq -i TR_TF.unmapped.sorted.bam -fq TR_TF.final.unmapped.R1.fq  -fq2 TR_TF.final.unmapped.R2.fq 2>/dev/null

    # 3) Extract mapped reads of the finalmapping (these are also, by definition, rotifer reads
    MAPPATH="$CLONE/finalmapping/"
    INFILE2="$MAPPATH$CLONE"_"$SAMPLENO"_to_VBCFpol.sorted.bam""
    samtools view -u -F4 $INFILE2 > TR_TF.mapped.bam # extract only reads that map to the reference genome 
    samtools sort -n TR_TF.mapped.bam -o TR_TF.mapped.sorted.bam
    bedtools bamtofastq -i TR_TF.mapped.sorted.bam -fq TR_TF.final.mapped.R1.fq  -fq2 TR_TF.final.mapped.R2.fq 2>/dev/null

    # 4) Merge these "rotifer reads" with mapped reads of the final assembly
    OUTREADS1="$OUTMAIN$CLONE"_"$SAMPLENO.finalreads.R1.fq"
    OUTREADS2="$OUTMAIN$CLONE"_"$SAMPLENO.finalreads.R2.fq"
    cat TR_TF.final.mapped.R1.fq TR_TF.final.unmapped.R1.fq > $OUTREADS1
    cat TR_TF.final.mapped.R2.fq TR_TF.final.unmapped.R2.fq > $OUTREADS2
    
    # Count reads and BPs after removal of tetraselmis
    FASTA_R1="TR_TF.R1.fasta"
    FASTA_R2="TR_TF.R2.fasta"
    COL[1]=$CLONE
    COL[2]=$SAMPLENO
    sed -n '1~4s/^@/>/p;2~4p' $OUTREADS1 > $FASTA_R1  # Convert to FASTA-format
    sed -n '1~4s/^@/>/p;2~4p' $OUTREADS2 > $FASTA_R2
    N_READS_R1=$(grep -c '^>' $FASTA_R1)  # number of reads in multi-fasta
    N_READS_R2=$(grep -c '^>' $FASTA_R2)
    (( COL[3] = $N_READS_R1 + $N_READS_R2 ))   # sum-up results of both mates
    TOTAL_BP_R1=$(grep -v '>' $FASTA_R1 | tr -d '\n' | wc -c) # total number of bases in fasta
    TOTAL_BP_R2=$(grep -v '>' $FASTA_R2 | tr -d '\n' | wc -c)
    (( COL[4] = $TOTAL_BP_R1 + $TOTAL_BP_R2 ))
    echo ${COL[*]} >> "SCRIPTS/after_tetra_removal.report.csv"

    # Run fastqc on finalreads2 
    FQCPATH="$OUTMAIN/fastqc"
	mkdir -p $FQCPATH
	fastqc --threads 8 $OUTREADS1 $OUTREADS2 -o $FQCPATH
	
    # Run jellyfish on finalreads2
    COMBIFQ="TR_TF.finalreads.combined.fq"
	cat $OUTREADS1 $OUTREADS2 > $COMBIFQ
	JELLYPATH="$OUTMAIN/jellyfish"
	mkdir -p $JELLYPATH	
	JELLYFILE="$JELLYPATH/$CLONE"_"$SAMPLENO.jf"
	JELLYHISTO="$JELLYPATH/$CLONE"_"$SAMPLENO.histo"
	jellyfish count -C -m 21 -s 1000000000 -t 8 $COMBIFQ -o $JELLYFILE
	jellyfish histo -t 1 $JELLYFILE > $JELLYHISTO
   
    # Clear all temporary files
    rm tffs1.bam tffs2.bam tffs3.bam TR_TF*
    
    # Compression of read-files at the end of the run
    echo "Compressing read files..."    
    gzip $OUTREADS1 $OUTREADS2 

done < SCRIPTS/datasets/OHJ_all.csv


IFS=$OLDIFS
