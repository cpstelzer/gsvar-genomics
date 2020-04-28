#!/bin/bash

# Running kraken2 on mapped reads (to determine the false discovery
# rate of kraken2 db with the protozoa database)

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

INMAIN="/aqtrim/mapped/"

cd ..  # escape "contaminants" folder

module load kraken2/20200325

while read CLONE SAMPLENO FRAGSIZE

do

	echo "Currently processing  $CLONE"	
    
    INGZ1="$CLONE$INMAIN$CLONE"_"$SAMPLENO.mapped.R1.fq.gz"
    INGZ2="$CLONE$INMAIN$CLONE"_"$SAMPLENO.mapped.R2.fq.gz"
    
    echo "Decompressing data..."
    gunzip $INGZ1 $INGZ2

    IN1="$CLONE$INMAIN$CLONE"_"$SAMPLENO.mapped.R1.fq"
    IN2="$CLONE$INMAIN$CLONE"_"$SAMPLENO.mapped.R2.fq"
    

    # Run kraken2 analysis on mapped reads to determine false positive rate of contaminant discovery    
    echo "Run kraken2 analysis on mapped reads (against standard database)..."    
    REPORT_MAPPED="$CLONE"_"$SAMPLENO.mapped.proto_db.log"
    kraken2 --threads 8 --db /blastdb/kraken2/protozoa --report $REPORT_MAPPED --paired $IN1 $IN2  >/dev/null
    mv $REPORT_MAPPED contaminants/kraken2/false_positives3/

    echo "Compressing data..."
    gzip $IN1 $IN2

done < SCRIPTS/datasets/OHJ_all.csv


IFS=$OLDIFS
