#!/bin/bash

# Script to collect fastqc-data on different steps during short-read preprocessing
# from automatically generated fastqc reports

# CP-Stelzer (Notebook, page 84, 29.4.2020)

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS



declare -a SUMPARS=("Basic Statistics"
	"Per base sequence quality"
	"Per tile sequence quality"
	"Per sequence quality scores"
	"Per base sequence content"
	"Per sequence GC content"
	"Per base N content"
	"Sequence Length Distribution"
	"Sequence Duplication Levels"
	"Overrepresented sequences"
	"Adapter Content"
	"Kmer Content")

#MAINPATH="rawdata"
#MAINPATH="aqtrim"
MAINPATH="finalreads"

OUTREPORT="SCRIPTS/fastqc.report.$MAINPATH.csv"

cd .. # escape SCRIPTS folder
mkdir -p TEMP_FQC # create a temporary folder to store unzipped fastqc reports

declare -a RID=("R1" "R2") # Counter for R1 or R2 files

while read CLONE SAMPLENO FRAGSIZE



do

    # Extracting fastqc files to TEMP_FQC folder (will be deleted after the run)
    
    # FILE1="$CLONE"/"$MAINPATH"/fastqc/"$CLONE"_"$SAMPLENO".R1_fastqc.zip""        # for rawreads-stage of preprocessing
    # FILE2="$CLONE"/"$MAINPATH"/fastqc/"$CLONE"_"$SAMPLENO".R2_fastqc.zip""
    # FILE1="$CLONE"/"$MAINPATH"/fastqc/"$CLONE"_"$SAMPLENO".R1.trimmed_fastqc.zip""  # for aqtrim-stage of preprocessing
    # FILE2="$CLONE"/"$MAINPATH"/fastqc/"$CLONE"_"$SAMPLENO".R2.trimmed_fastqc.zip""
    FILE1="$CLONE"/"$MAINPATH"/fastqc/"$CLONE"_"$SAMPLENO".roti-mito.dedup.R1_fastqc.zip""  # for finalreads-stage of preprocessing
    FILE2="$CLONE"/"$MAINPATH"/fastqc/"$CLONE"_"$SAMPLENO".roti-mito.dedup.R2_fastqc.zip""
    


    unzip $FILE1 -d TEMP_FQC/
    unzip $FILE2 -d TEMP_FQC/
    #FQCS1=""TEMP_FQC/"$CLONE"_"$SAMPLENO".R1_fastqc"" # new path of extracted folder (rawdata)
    #FQCS2=""TEMP_FQC/"$CLONE"_"$SAMPLENO".R2_fastqc""
    #FQCS1=""TEMP_FQC/"$CLONE"_"$SAMPLENO".R1.clean_fastqc"" # new path of extracted folder (aqtrim)
    #FQCS2=""TEMP_FQC/"$CLONE"_"$SAMPLENO".R2.clean_fastqc""
    FQCS1=""TEMP_FQC/"$CLONE"_"$SAMPLENO".roti-mito.dedup.R1_fastqc"" # new path of extracted folder (finalreads)
    FQCS2=""TEMP_FQC/"$CLONE"_"$SAMPLENO".roti-mito.dedup.R2_fastqc""



    for i in "${RID[@]}"
    do
        
        CURRENT="$i"
    
        # Switching between R1 and R2
        if [ "$CURRENT" = "R1" ]; then
            echo "Currently reading from R1"
            FILENAME="$CLONE"_"$SAMPLENO"."$MAINPATH"."$CURRENT"       
            SUMMARYFILE="$FQCS1"/summary.txt""
            DATAFILE="$FQCS1"/fastqc_data.txt"" 
        else
            echo "Currently reading from R2"
            FILENAME="$CLONE"_"$SAMPLENO"."$MAINPATH"."$CURRENT"       
            SUMMARYFILE="$FQCS2"/summary.txt""
            DATAFILE="$FQCS2"/fastqc_data.txt"" 
        fi

	    
	    COL[1]=$FILENAME
	    COL[2]=$(grep -Po '[\w]{4}(?=\tBasic Statistics)' $SUMMARYFILE)
	    COL[3]=$(grep -Po '[\w]{4}(?=\tPer base sequence quality)' $SUMMARYFILE)
	    COL[4]=$(grep -Po '[\w]{4}(?=\tPer tile sequence quality)' $SUMMARYFILE)
	    COL[5]=$(grep -Po '[\w]{4}(?=\tPer sequence quality scores)' $SUMMARYFILE)
	    COL[6]=$(grep -Po '[\w]{4}(?=\tPer base sequence content)' $SUMMARYFILE)
	    COL[7]=$(grep -Po '[\w]{4}(?=\tPer sequence GC content)' $SUMMARYFILE)
	    COL[8]=$(grep -Po '[\w]{4}(?=\tPer base N content)' $SUMMARYFILE)
	    COL[9]=$(grep -Po '[\w]{4}(?=\tSequence Length Distribution)' $SUMMARYFILE)
	    COL[10]=$(grep -Po '[\w]{4}(?=\tSequence Duplication Levels)' $SUMMARYFILE)
	    COL[11]=$(grep -Po '[\w]{4}(?=\tOverrepresented sequences)' $SUMMARYFILE)
	    COL[12]=$(grep -Po '[\w]{4}(?=\tAdapter Content)' $SUMMARYFILE)
	    COL[13]=$(grep -Po '[\w]{4}(?=\tKmer Content)' $SUMMARYFILE)
	    COL[14]=$(grep -Po '(?<=Total Sequences\s)[\d]*' $DATAFILE)
	    COL[15]=$(grep -Po '(?<=Sequence length\s)[\d|\-]*' $DATAFILE)
	    COL[16]=$(grep -Po '(?<=%GC\s)[\d]*' $DATAFILE)
	    echo ${COL[*]}
	    echo ${COL[*]} >> $OUTREPORT

    done # for loop

done < SCRIPTS/datasets/OHJ_all.csv

rm -r TEMP_FQC # delete the temporary folder

IFS=$OLDIFS



