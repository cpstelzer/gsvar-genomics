#!/bin/bash

# Map (mapped or unmapped) finalreads to the genome of Tetraselmis striata, the closest relative of 
# our food alga with a published genome
# see Notebook page 91 and ff. (13.5.2020)

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

#cd tetraselmis
#bowtie2-build GCA_006384855.1_TSEL_PacBio_SMRT_genomic.fna  tetra_striata_IDX # Building index (needs to be done only once)
#cd ..
TETRAIDX="/contaminants/tetraselmis/tetra_striata_IDX" # Path to the index of the Tetraselmis striata genome

cd ..  # escape contaminants folder


module load bedtools/2.27.1

while read CLONE SAMPLENO FRAGSIZE

do
    echo "Currently processing  $CLONE"	
	
    MAPPATH="$CLONE/finalmapping/"

    mkdir -p "$MAPPATH/tetraselmis/"

    INFILE="$MAPPATH$CLONE"_"$SAMPLENO"_to_VBCFpol.sorted.bam""

    # Filtering out reads that did NOT map in the final alignment...
    samtools view -u -f 4 -F 264 $INFILE > tmps1.bam       # An unmapped read whose mate is mapped
    samtools view -u -f 8 -F 260 $INFILE > tmps2.bam       # A mapped read who’s mate is unmapped
    samtools view -u -f 12 -F 256 $INFILE > tmps3.bam      # Both reads of the pair are unmapped
    samtools merge temp.unmapped.bam tmps1.bam tmps2.bam tmps3.bam
    samtools sort -n temp.unmapped.bam -o temp.unmapped.sorted.bam
    bedtools bamtofastq -i temp.unmapped.sorted.bam -fq temp.final.unmapped.R1.fq  -fq2 temp.final.unmapped.R2.fq 2>/dev/null
    
    # ...map those reads to the Tetraselmis genome
    OUTLOG1="$MAPPATH"/tetraselmis/"$CLONE"_"$SAMPLENO".unmapped_finalreads_2_tetra.output.log""
    OUTBAM="$MAPPATH"/tetraselmis/"$CLONE"_"$SAMPLENO".unmapped_finalreads_2_tetra.sorted.bam""
    (bowtie2 -p 8 -x $TETRAIDX -X $FRAGSIZE -1 temp.final.unmapped.R1.fq -2 temp.final.unmapped.R2.fq -S temp.unmapped2tetra.sam)2> $OUTLOG1
	samtools view -S -h -@ 8 temp.unmapped2tetra.sam | samtools sort -o $OUTBAM  # for later use, if contamination with tetra is real
    # Unmapped final reads that map to Tetraselmis --> potential contaminants
    READS1=$(grep -Po '(?<=\s)[0-9]*(?=\s[\(|p|m])' $OUTLOG1)
    COL[1]=$CLONE
    COL[2]=$SAMPLENO    
    COL[3]=$READS1
    echo ${COL[*]} >> contaminants/map2tetra.potential_contam.report.csv

    # Filtering out reads that DID map in the final alignment
    samtools view -u -F4 $INFILE > temp.mapped.bam # extract only reads that map to the reference genome 
    samtools sort -n temp.mapped.bam -o temp.mapped.sorted.bam
    bedtools bamtofastq -i temp.mapped.sorted.bam -fq temp.final.mapped.R1.fq  -fq2 temp.final.mapped.R2.fq 2>/dev/null

    # ...map those reads to the Tetraselmis genome
    OUTLOG2="$MAPPATH"/tetraselmis/"$CLONE"_"$SAMPLENO".mapped_finalreads_2_tetra.output.log""
    (bowtie2 -p 8 -x $TETRAIDX -X $FRAGSIZE -1 temp.final.mapped.R1.fq -2 temp.final.mapped.R2.fq -S temp.mapped2tetra.sam)2> $OUTLOG2
	# samtools view -S -h -@ 8 temp.mapped2tetra.sam | samtools sort -o temp.mapped2tetra.sorted.bam     # this is actually not needed!
    # Mapped final reads that map to Tetraselmis --> false positives
    READS2=$(grep -Po '(?<=\s)[0-9]*(?=\s[\(|p|m])' $OUTLOG2)
    CCOL[1]=$CLONE
    CCOL[2]=$SAMPLENO    
    CCOL[3]=$READS2
    echo ${CCOL[*]} >> contaminants/map2tetra.false_positives.report.csv

    # Clear all temporary files
    rm tmps1.bam tmps2.bam tmps3.bam temp*

done < SCRIPTS/datasets/OHJ_all.csv


IFS=$OLDIFS

