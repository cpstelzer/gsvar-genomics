#!/bin/bash

# Decontamination pipeline according to the schematic in notebook page 56 (30.3.2020) 
# Data dependencies: Requires mappings to Pseudomonas toyotomiensis already done (map2pseudo.sh)
# and BAM files of these alignments moved to subfolder "map2pseudoT" within $CLONE/aqtrim/unmapped.

# Alternatively, one could skip the Pseudomonas-step and directly use the unmapped reads from the
# original BAM file...

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

INMAIN="/aqtrim/unmapped/map2pseudoT/"
K2OUT="/aqtrim/unmapped/kraken2/"
ROTIREADS="/rotireads/"
ASSEMBLY="PseudoT"
KRAKEN2_STD="/blastdb/kraken2/standard-no-human"  # Rich made an extra DB without human genome

cd ..  # escape "contaminants" folder
module load bedtools/2.27.1
module load kraken2/20200325

while read CLONE SAMPLENO FRAGSIZE

do

	echo "Currently processing  $CLONE"	
	
	# 1) Extract reads that don't map to Pseudomonas toyotomiensis 
    INFILE="$CLONE$INMAIN$CLONE"_"$SAMPLENO"_to_"$ASSEMBLY".sorted.bam""
    # Filtering out reads that did not map (code from novocraft.com)
    samtools view -u -f 4 -F 264 $INFILE > tmps1.bam       # An unmapped read whose mate is mapped
    samtools view -u -f 8 -F 260 $INFILE > tmps2.bam       # A mapped read who’s mate is unmapped
    samtools view -u -f 12 -F 256 $INFILE > tmps3.bam      # Both reads of the pair are unmapped
    samtools merge unmapped.bam tmps1.bam tmps2.bam tmps3.bam
    samtools sort -n unmapped.bam -o unmapped.sorted.bam
	bedtools bamtofastq -i unmapped.sorted.bam -fq unmapped.R1.fq  -fq2 unmapped.R2.fq

    # 2) Run kraken2 analysis with the standard database (i.e., libraries of bacteria, archaea, viruses, UniVec_core)
    echo "Run kraken2 analysis with the standard database..."  
    K2FOLDER="$CLONE$K2OUT"
    mkdir -p $K2FOLDER 
    REPORT_UNMAPPED_STD="$CLONE"_"$SAMPLENO.kraken2.std_db.log"
    kraken2 --threads 10 --db $KRAKEN2_STD --report $REPORT_UNMAPPED_STD \
    --paired --classified-out cseqs#.fq --unclassified-out useqs#.fq unmapped.R1.fq unmapped.R2.fq >/dev/null 
    
    # Rename and store classified reads
    CLASS_STD_R1="$CLONE"_"$SAMPLENO.kraken2.std_db.R1.fq"
    CLASS_STD_R2="$CLONE"_"$SAMPLENO.kraken2.std_db.R2.fq"
    mv cseqs_1.fq $CLASS_STD_R1
    mv cseqs_2.fq $CLASS_STD_R2
    mv $CLASS_STD_R1 $K2FOLDER 
    mv $CLASS_STD_R2 $K2FOLDER 
    mv $REPORT_UNMAPPED_STD $K2FOLDER
    
    mv useqs_1.fq remaining.R1.fq # rename unclassified reads for step 3
    mv useqs_2.fq remaining.R2.fq
    
    
    # 3)  Run kraken2 analysis with the protozoa database
    echo "Run kraken2 analysis with the protozoa database..."
    REPORT_UNMAPPED_PROTO="$CLONE"_"$SAMPLENO.kraken2.proto_db.log"  
    kraken2 --threads 10 --db /blastdb/kraken2/protozoa --report $REPORT_UNMAPPED_PROTO \
        --paired --classified-out cseqs#.fq --unclassified-out useqs#.fq remaining.R1.fq remaining.R2.fq >/dev/null
    
    # Rename and store the classified reads    
    CLASS_STD_R1="$CLONE"_"$SAMPLENO.kraken2.proto_db.R1.fq"
    CLASS_STD_R2="$CLONE"_"$SAMPLENO.kraken2.proto_db.R2.fq"
    mv cseqs_1.fq $CLASS_STD_R1
    mv cseqs_2.fq $CLASS_STD_R2
    mv $CLASS_STD_R1 $K2FOLDER   # save classified reads 
    mv $CLASS_STD_R2 $K2FOLDER
    mv $REPORT_UNMAPPED_PROTO $K2FOLDER

    # 4) Declare all the remaining, still unclassified reads as "unmapped rotifer reads"
    # Rename and store these reads & combine them with the mapped rotifer reads
    # to get the "total rotifer reads"
    mkdir -p "$CLONE$ROTIREADS"

    # Retrieve mapped reads from the original BAM file
    # echo "Retrieve mapped reads from the original BAM file..."
    # BAMFILE="$CLONE"/mapping/"$CLONE"_"$SAMPLENO"_to_VBCFpol.sorted.bam"" # the original alignment to the reference
    # samtools view -u -F4 $BAMFILE > mapped.bam # extract only reads that map to the reference genome 
    # samtools sort -n mapped.bam -o mapped.sorted.bam
	# bedtools bamtofastq -i mapped.sorted.bam -fq mapped.R1.fq  -fq2 mapped.R2.fq 2>/dev/null

    # Load mapped reads directly from stored fq.gz file (alternative to the above extraction from BAM file)
    MR1="$CLONE"/aqtrim/mapped/"$CLONE"_"$SAMPLENO.mapped.R1.fq.gz"
    MR2="$CLONE"/aqtrim/mapped/"$CLONE"_"$SAMPLENO.mapped.R2.fq.gz"
    zcat $MR1 > mapped.R1.fq
    zcat $MR2 > mapped.R2.fq
        
    # Combine mapped reads and "unmapped rotifer reads"
    ROTI_R1="$CLONE$ROTIREADS$CLONE"_"$SAMPLENO.rotireads.R1.fq"
    ROTI_R2="$CLONE$ROTIREADS$CLONE"_"$SAMPLENO.rotireads.R2.fq"
    cat mapped.R1.fq useqs_1.fq > $ROTI_R1
    cat mapped.R2.fq useqs_2.fq > $ROTI_R2
    gzip $ROTI_R1
    gzip $ROTI_R2

    # Run kraken2 analysis on mapped reads to determine false positive rate of contaminant discovery    
    echo "Run kraken2 analysis on mapped reads (against standard database)..."    
    REPORT_MAPPED="$CLONE"_"$SAMPLENO.mapped.std_db.log"
    kraken2 --threads 10 --db $KRAKEN2_STD --report $REPORT_MAPPED --paired mapped.R1.fq mapped.R2.fq  >/dev/null
    mv $REPORT_MAPPED contaminants/kraken2/false_positives3/

    # Additionally store "unmapped rotifer reads"
    UNCLASS_R1="$CLONE"_"$SAMPLENO.kraken2.unclass.R1.fq"
    UNCLASS_R2="$CLONE"_"$SAMPLENO.kraken2.unclass.R2.fq"
    mv useqs_1.fq $UNCLASS_R1
    mv useqs_2.fq $UNCLASS_R2
    mv $UNCLASS_R1 $K2FOLDER
    mv $UNCLASS_R2 $K2FOLDER

    # Clear all temporary files
    rm tmps*.bam unmapped.bam unmapped.sorted.bam \
       unmapped.R1.fq unmapped.R2.fq mapped.R1.fq mapped.R2.fq \
       remaining.R1.fq remaining.R2.fq \
       # mapped.bam mapped.sorted.bam 
   
    
done < SCRIPTS/datasets/OHJ_all.csv


IFS=$OLDIFS
