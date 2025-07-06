#!/bin/bash

set -xvu

##### Script to map met1 HiFi reads to WT genome
# USAGE: #	
# ./Scripts/script_aln-HiFi.sh \
#  ./ALN \
#  ${SP} \
#  ./Tsu-0.reference.fa
#  32

# for SP in `cat ./list.SINGLE-PLANTS.txt` ; do
#  ./Scripts/script_aln-HiFi.sh \
#   ./ALN \
#   ${SP} \
#   ./Tsu-0.reference.fa
#   32
# done

# Parse parameteres
WORKDIR=$1
SINGLEPLANT=$2
REF=$3
CORES=${4:-${NSLOTS:-12}}


# Declare other variables

# seq file ./Reads/${SINGLEPLANT}.dc_q20.fastq.gz
# reseq file ./Reads/${SINGLEPLANT}r.dc_q20.fastq.gz
reads=./Data/Reads/${SINGLEPLANT}merged.dc_q20.fastq.gz

mkdir -p $WORKDIR/merylDB

# Virtual env activation:
# containing winnowmap v2.03; samtools v1.16.1
## source activate /ebio/abt6_projects/met1_somatic_transpositions/conda/so_ALNtransposition

# Start
echo ""
echo "Start"
date

# merging reads from sequencing and resequencing SMRTcells
## WT comes as a single SMRTcell
zcat ./Data/Reads/${SINGLEPLANT}.dc_q20.fastq.gz ./Data/Reads/${SINGLEPLANT}r.dc_q20.fastq.gz > $reads

# winnowmap mapping reads to WT genome
meryl count k=15 output $WORKDIR/merylDB $REF
meryl print greater-than distinct=0.9998 $WORKDIR/merylDB > $WORKDIR/merylDB/repetitive_k15.txt

winnowmap -W repetitive_k15.txt -Y -L  -t ${CORES} -ax map-pb ${REF} ${reads} | samtools view -h -b -  > $WORKDIR/$SINGLEPLANT-aln2-WT-Chr_winnowmap.bam  
samtools sort -O bam -o $WORKDIR/$SINGLEPLANT-aln2-WT-Chr_winnowmap.sorted.bam $WORKDIR/$SINGLEPLANT-aln2-WT-Chr_winnowmap.bam


samtools addreplacerg -w -r "@RG\tID:${SINGLEPLANT}\tLB:${SINGLEPLANT}\tPL:sequelII_hifi_v3\tSM:${SINGLEPLANT}" \
		-@ $CORES -o $WORKDIR/$SINGLEPLANT-aln2-WT-Chr_winnowmap.sortedRG.bam \
			$WORKDIR/$SINGLEPLANT-aln2-WT-Chr_winnowmap.sorted.bam -O BAM &&
			samtools index $WORKDIR/$SINGLEPLANT-aln2-WT-Chr_winnowmap.sortedRG.bam

# # Calculating coverage
# covoutput=$WORKDIR/coverage
# mkdir -p  $covoutput

# ## making 500bp windows
# bioawk -c fastx '{print $name, length($seq)}' $REF > $covoutput/4makingwindows.txt
# bedtools makewindows -g $covoutput/4makingwindows.txt -w 500 > $covoutput/w500.bed

# ## Calculating coverage at 500bp windows
# mosdepth --fast-mode --by $covoutput/w500.bed $covoutput/$SINGLEPLANT.w500 $WORKDIR/$SINGLEPLANT-aln2-WT_winnowmap.sortedRG.bam

# ## 100 windows with higher coverage
# zcat $covoutput/$SINGLEPLANT.w500.regions.bed.gz | sort -g -k4 -r | head -100 > $covoutput/highest-100_cov.txt


# rm intermediary files
rm $WORKDIR/$SINGLEPLANT-aln2-WT_winnowmap.sorted.bam 

# End
date
echo "End"
echo ""

# Virtual env deactivation
## conda deactivate