#!/bin/bash

OLDIFS=$IFS     # save the existing field separator
IFS=";"

cd ..  # escape SCRIPTS folder
module load bbtools/38.34
INMAIN="/rawdata/"
OUTMAIN="/aqtrim/"
OUTMAPP="/mapping/"
ASSEMBLY="VBCFpol"   # Polished genome. Alternatively, VBCF (=unpolished version)

while read CLONE SAMPLENO FRAGSIZE

do
	echo "Currently processing  $CLONE"	

	mkdir -p "$CLONE$OUTMAIN"
	mkdir -p "$CLONE$OUTMAPP"	

	INFILE1="$CLONE$INMAIN$CLONE"_"$SAMPLENO.R1.fq"
	INFILE2="$CLONE$INMAIN$CLONE"_"$SAMPLENO.R2.fq"
	# OUTFILE1="$CLONE$OUTMAIN$CLONE"_"$SAMPLENO.R1.trimmed.fq" # OLD version --> decactivated
	# OUTFILE2="$CLONE$OUTMAIN$CLONE"_"$SAMPLENO.R2.trimmed.fq"
	OUTFILE1="$CLONE$OUTMAIN$CLONE"_"$SAMPLENO.trimmed.R1.fq" # new version 3.4.2020 !!!
        OUTFILE2="$CLONE$OUTMAIN$CLONE"_"$SAMPLENO.trimmed.R2.fq"

	#echo "$INFILE1"
	#echo "$INFILE2"
	#echo "$OUTFILE1"
	#echo "$OUTFILE2"
	
	echo "Run bbduk..."
	OUTREPORT="$CLONE$OUTMAIN$CLONE"_"$SAMPLENO.report.txt"
	bbduk.sh in=$INFILE1 in2=$INFILE2 out=$OUTFILE1 out2=$OUTFILE2 \
	ref=adaptors.fa k=23 ktrim=n mink=11 hdist=1 tpe tbo  qtrim=rl trimq=20 maq=10 minlen=40 &> $OUTREPORT

        FQCPATH="$CLONE/aqtrim/fastqc"
	mkdir -p $FQCPATH
	fastqc --threads 8 $OUTFILE1 $OUTFILE2 -o $FQCPATH
	echo $FQCPATH

	echo "Run jellyfish"
	COMBIFQ="$CLONE$OUTMAIN$CLONE"_"$SAMPLENO.all.trimmed.fq"
	cat $OUTFILE1 $OUTFILE2 > $COMBIFQ
	JELLYPATH="$CLONE/aqtrim/jellyfish"
	mkdir -p $JELLYPATH	
	JELLYFILE="$JELLYPATH/$CLONE"_"$SAMPLENO.jf"
	JELLYHISTO="$JELLYPATH/$CLONE"_"$SAMPLENO.histo"
	jellyfish count -C -m 21 -s 1000000000 -t 10 $COMBIFQ -o $JELLYFILE
	jellyfish histo -t 1 $JELLYFILE > $JELLYHISTO

	
	echo "Run Bowtie2 mapping..."
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
	
	(bowtie2 -p 8 -x $GINDEX -X $FRAGSIZE -1 $OUTFILE1 -2 $OUTFILE2 -S $OUTSAM)2> $OUTLOG
	samtools view -S -h -@ 10 $OUTSAM | samtools sort -o $OUTBAM
	MAPPATH="$CLONE$OUTMAPP"
	mkdir $MAPPATH
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
	rm $COMBIFQ
	gzip $INFILE1
	gzip $INFILE2
	gzip $OUTFILE1
	gzip $OUTFILE2

	#echo $OUTSAM
	#echo $OUTLOG
	#echo $OUTBAM
	#echo $MAPPATH
	#echo $BAMFILE

	#echo $COMBIFQ
	#echo $JELLYPATH
	#echo $JELLYFILE
	#echo $JELLYHISTO
	

done < SCRIPTS/rerun.csv


IFS=$OLDIFS
