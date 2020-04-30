#!/bin/bash

# script for collecting information on mapping of Illumina reads to Pacbio assembly
# see page 36 in notebook (9.3.2020) for more details

ASSEMBLY="VBCFpol"
INMAIN="finalmapping"
OUTMAIN="SCRIPTS/$ASSEMBLY.$INMAIN.report.csv"

# Information on Headers
# total: 	total no of paired reads
# ac0t:		pairs that aligned concordantly 0 times
# ac1t		pairs that aligned concordantly exactly 1 time (i.e., exactly 1 place in the assembly)
# acg1t		pairs that aligned concordantly more than 1 time
# ac0t		see above
# ad1t		pairs that aligned discordantly exactly 1 time
# acd0t		pairs that did not align concordantly nor discordantly
# mmp		mates that make up the acd0t
# sr0t		singletons that do not align anywhere
# sr1t		singletons that align exactly 1 time (but the other mate does not align)
# srg1t		singletons that align more than 1 time	(mut the other mate does not align)


cd ..

OLDIFS=$IFS
NEWIFS=";"
IFS=$NEWIFS

while read CLONE VBCFNO ADDCOL

do
     INFILE="$CLONE/$INMAIN/$CLONE"_"$VBCFNO"_to_"$ASSEMBLY".output.log""
     READS=$(grep -Po '(?<=\s)[0-9]*(?=\s[\(|p|m])' $INFILE)
     COL[1]=$CLONE
     COL[2]=$VBCFNO
     COL[3]=$READS
     echo ${COL[*]} 
     IFS=$OLDIFS
     echo ${COL[*]} >> $OUTMAIN
     IFS=$NEWIFS

done < SCRIPTS/datasets/OHJ_all.csv

IFS=$OLDIFS
