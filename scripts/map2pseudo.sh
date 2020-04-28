#!/bin/bash

# Mapping previously unmapped reads against Pseudomonas toyotomiensis genome
# see Notebook page 45 (19.3.2020)
# see also script move2pseudoFolder.sh (page 56, 31.3.2020) which was run afterwards


OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

cd ..  # escape "contaminants" folder
INMAIN="/aqtrim/unmapped/"
GINDEX="contaminants/quast_results/Pseudo_toy_IDX"   # Index of contaminant genome Pseudomonas 
ASSEMBLY="PseudoT"

while read CLONE SAMPLENO FRAGSIZE

do
	echo "Currently processing  $CLONE"	

	INFILE1="$CLONE$INMAIN$CLONE"_"$SAMPLENO.unmapped.R1.fq.gz"
	INFILE2="$CLONE$INMAIN$CLONE"_"$SAMPLENO.unmapped.R2.fq.gz"

	echo "Run Bowtie2 mapping..."
	OUTSAM="$CLONE$INMAIN$CLONE"_"$SAMPLENO"_to_"$ASSEMBLY".sam""
	OUTLOG="$CLONE$INMAIN$CLONE"_"$SAMPLENO"_to_"$ASSEMBLY".output.log""
	OUTBAM="$CLONE$INMAIN$CLONE"_"$SAMPLENO"_to_"$ASSEMBLY".sorted.bam""
	
	(bowtie2 -p 8 -x $GINDEX -X $FRAGSIZE -1 $INFILE1 -2 $INFILE2 -S $OUTSAM)2> $OUTLOG
	samtools view -S -h -@ 10 $OUTSAM | samtools sort -o $OUTBAM
	MAPPATH="$CLONE$INMAIN"
	BAMFILE="$CLONE"_"$SAMPLENO"_to_"$ASSEMBLY".sorted.bam""
	cd $MAPPATH
	samtools index $BAMFILE
	cd ../../..

	# Extract mapping information
	OUTREPORT="contaminants/$ASSEMBLY.mapping.report.csv"
    # Information on Headers
    # total:        total no of paired reads
    # ac0t:         pairs that aligned concordantly 0 times
    # ac1t          pairs that aligned concordantly exactly 1 time (i.e., exactly 1 place in the assembly)
    # acg1t         pairs that aligned concordantly more than 1 time
    # ac0t          see above
    # ad1t          pairs that aligned discordantly exactly 1 time
    # acd0t         pairs that did not align concordantly nor discordantly
    # mmp           mates that make up the acd0t
    # sr0t          singletons that do not align anywhere
    # sr1t          singletons that align exactly 1 time (but the other mate does not align)
    # srg1t         singletons that align more than 1 time  (mut the other mate does not align)
    READS=$(grep -Po '(?<=\s)[0-9]*(?=\s[\(|p|m])' $OUTLOG)
    COL[1]=$CLONE
    COL[2]=$SAMPLENO
    COL[3]=$READS
    echo ${COL[*]}
    IFS=$OLDIFS
    echo ${COL[*]} >> $OUTREPORT
    IFS=$NEWIFS

	# Extract unmapped reads that do NOT align to the contaminant genome
    
    # Delete and zip intermediate files to save space
	rm $OUTSAM
    #rm $OUTBAM
	
	

done < SCRIPTS/datasets/OHJ7_mbl2019.csv


IFS=$OLDIFS
