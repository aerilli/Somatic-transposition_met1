#!/bin/bash


# amovilli on 2023-10-09

# Script to remove TEs overlapping with other repeats (centromere, telomere, rDNA, organellar)

## EXAMPLE ## /ebio/abt6_projects/met1_somatic_transpositions/code/script_RemovingOtherRepeatsOverlaps.sh \
##              /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/5_TEAnno/EDTA \
##              /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/5_TEAnno/EDTA/Tsu-0.scaffolds_contigs.v2.fa.mod.EDTA.intact.gff3


set -xv

# Parse parameteres
WORKDIR=$1
INPUT=$2

# Declare other variables
output=$WORKDIR/`basename $INPUT ".gff3"`.no-other-repeats.gff3
otherRepeats="/ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/3_annotation_REPEATS/Chr_scaffolds/Tsu-0/Tsu-0.Repeats_merged.gff"


# Start
echo ""
echo "Start"
date

# TEs have to overlap at least 50% to be removed
bedtools intersect -header -v -f 0.5 -a $INPUT -b $otherRepeats > $output


# End
echo ""
echo "Start"
date