#!/bin/bash

# Mapping the finally prepared reads ($CLONE"_"$SAMPLENO.roti-mito.dedup.R[1|2].fq to the reference genome (VBCFpol)
# See Notebook page 78, 24.4.2020

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

cd ..  # escape SCRIPTS folder

INMAIN="/finalreads/"
OUTMAPP="/finalmapping/"
ASSEMBLY="VBCFpol"   # Polished genome. Alternatively, VBCF (=unpolished version)

while read CLONE SAMPLENO FRAGSIZE

do
	echo "Currently processing input-files of  $CLONE"_"$SAMPLENO"	

	INGZ1="$CLONE$INMAIN$CLONE"_"$SAMPLENO.roti-mito.dedup.R1.fq.gz"
    INGZ2="$CLONE$INMAIN$CLONE"_"$SAMPLENO.roti-mito.dedup.R2.fq.gz"

    gunzip $INGZ1 $INGZ2

    INFILE1="$CLONE$INMAIN$CLONE"_"$SAMPLENO.roti-mito.dedup.R1.fq"
	INFILE2="$CLONE$INMAIN$CLONE"_"$SAMPLENO.roti-mito.dedup.R2.fq"
	
	echo "Run Bowtie2 mapping..."

    MAPPATH="$CLONE$OUTMAPP"
	mkdir -p $MAPPATH

	OUTSAM="$CLONE$OUTMAPP$CLONE"_"$SAMPLENO"_to_"$ASSEMBLY".sam""
	OUTLOG="$CLONE$OUTMAPP$CLONE"_"$SAMPLENO"_to_"$ASSEMBLY".output.log""
	OUTBAM="$CLONE$OUTMAPP$CLONE"_"$SAMPLENO"_to_"$ASSEMBLY".sorted.bam""
	
	if [ "$ASSEMBLY" == "VBCFpol" ]; then
    	   echo "...against the polished genome."
	   GINDEX="Baspl_polished_IDX"
	else
    	   echo "...against the unpolished genome."
	   GINDEX="Baspl_IDX"
	fi
	
	(bowtie2 -p 8 -x $GINDEX -X $FRAGSIZE -1 $INFILE1 -2 $INFILE2 -S $OUTSAM)2> $OUTLOG
	samtools view -S -h -@ 10 $OUTSAM | samtools sort -o $OUTBAM
	
    BAMFILE="$CLONE"_"$SAMPLENO"_to_"$ASSEMBLY".sorted.bam""
	cd $MAPPATH
	samtools index $BAMFILE
	cd ../..
	
	# Extract mapping information
	OUTCOV="$CLONE$OUTMAPP$CLONE"_"$SAMPLENO"_to_"$ASSEMBLY".50kbp.cov""
	# echo $OUTCOV
	samtools bedcov OHJ7i3n10pol_50kbpwin  $OUTBAM > $OUTCOV

    # Delete and zip intermediate files to save space
	rm $OUTSAM
	
    echo "Re-Compressing input files..."

	gzip $INFILE1
	gzip $INFILE2
		

done < SCRIPTS/datasets/OHJ_all_subset4.csv


IFS=$OLDIFS
