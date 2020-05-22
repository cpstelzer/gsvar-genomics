
#!/bin/bash

# Pipeline zur Extraktion der "primary absolute values" (siehe Notizbuch S. 62 & 63, 3.4.2020)

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

cd ..  # escape "contaminants" folder
module load bedtools/2.27.1

while read CLONE SAMPLENO FRAGSIZE

do # loop through clones/libraries

    echo "Currently processing  $CLONE"_"$SAMPLENO"	

    NCOUNT=0
    COL[1]=$CLONE
    COL[2]=$SAMPLENO
    
    while read VARNAME PATH2FILE FILEXT ZIPPED

    do  # loop through measured variables
        # VARNAME="total_reads"    
        # PATH2FILE="/aqtrim/"
        # FILEXT="trimmed"
        # ZIPPED=1
        NCOUNT=$NCOUNT+4
        
        FILE_GZ_R1="$CLONE$PATH2FILE$CLONE"_"$SAMPLENO"."$FILEXT".R1.fq.gz"" # compressed fq
        FILE_GZ_R2="$CLONE$PATH2FILE$CLONE"_"$SAMPLENO"."$FILEXT".R2.fq.gz""
        FILE_R1="$CLONE$PATH2FILE$CLONE"_"$SAMPLENO"."$FILEXT".R1.fq""       # uncompressed fq
        FILE_R2="$CLONE$PATH2FILE$CLONE"_"$SAMPLENO"."$FILEXT".R2.fq""
        FASTA_R1="$CLONE$PATH2FILE$CLONE"_"$SAMPLENO".R1.tmp""      # as (TEMPORARY) fasta-file
        FASTA_R2="$CLONE$PATH2FILE$CLONE"_"$SAMPLENO".R2.tmp""
        
        if [ $ZIPPED = "true" ]; then
            gunzip $FILE_GZ_R1 $FILE_GZ_R2
        fi
        
        # Convert FASTQ to FASTA
        sed -n '1~4s/^@/>/p;2~4p' $FILE_R1 > $FASTA_R1  # from: https://github.com/stephenturner/oneliners
        sed -n '1~4s/^@/>/p;2~4p' $FILE_R2 > $FASTA_R2

        # Extract no. of reads and totalBP
        N_READS_R1=$(grep -c '^>' $FASTA_R1)  # number of reads in multi-fasta
        N_READS_R2=$(grep -c '^>' $FASTA_R2)
        TOTAL_BP_R1=$(grep -v '>' $FASTA_R1 | tr -d '\n' | wc -c) # total number of bases in fasta
        TOTAL_BP_R2=$(grep -v '>' $FASTA_R2 | tr -d '\n' | wc -c)
        (( COL[$NCOUNT-1] = $N_READS_R1 + $N_READS_R2 ))   # sum-up results of both mates
        (( COL[$NCOUNT] = $TOTAL_BP_R1 + $TOTAL_BP_R2 ))

        # Clean up after data extraction
        gzip $FILE_R1 $FILE_R2; rm $FASTA_R1 $FASTA_R2
        
    done < contaminants/decon_vars.csv

    echo ${COL[*]}
    echo ${COL[*]} >> "contaminants/deconreport.csv" 


done < SCRIPTS/datasets/OHJ_all.csv


IFS=$OLDIFS
