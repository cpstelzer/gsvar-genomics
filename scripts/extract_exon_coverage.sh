#!/bin/bash

# Extracting coverage from regions annotated as "exon" in the polished B.asplanchnoidis assembly
# See notebook page 114, 5.6.2020

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

cd ..  # escape SCRIPTS folder
module load bbtools/38.34
INMAIN="/finalmapping/"
ASSEMBLY="VBCFpol"   # Polished genome. Alternatively, VBCF (=unpolished version)

while read CLONE SAMPLENO FRAGSIZE

do
	echo "Currently processing  $CLONE $SAMPLENO"	

    INBAM="$CLONE$INMAIN$CLONE"_"$SAMPLENO"_to_"$ASSEMBLY".sorted.bam""
	OUTCOV="$CLONE$INMAIN$CLONE"_"$SAMPLENO"_to_"$ASSEMBLY".exon.cov""

    IFS=$OLDIFS
	samtools bedcov VBCFpol_exon_win $INBAM > $OUTCOV
    # Note: brachionus.v2-round3.exwin.bed was renamed to 'OHJ7i3n10pol_exonwin'
    IFS=$NEWIFS

done < SCRIPTS/datasets/OHJ_all.csv


IFS=$OLDIFS

