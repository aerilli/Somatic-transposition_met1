#!/bin/bash

## I assume that  TEs of the same kind that are neighbouring each other are due to EDTA-derived fragmentation

INPUT=$1 # /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/5_TEanno/EDTA/Tsu-0.scaffolds_contigs.v2.fa.mod.EDTA.TEanno.no-other-repeats.gff3

output="${INPUT%.gff*}.no-parent.neighbours-merged.gff3"

echo "Start"
date
echo ""

for i in `cut -f9 ${INPUT} | sed 's/.*Name=//'| sed 's/;Classification.*//' | sed 's/AT[1-5]*TE[0-9]*_//' | sort | grep -v "^#" | uniq`; do
    grep -w "${i}" $INPUT | grep -v Parent|
    bedtools merge -i stdin -s -d 1 -c 2,3,9,6,7,8 -o distinct,distinct,distinct,distinct,distinct,distinct |
    awk '{OFS="\t"}{print $1,$4,$5,$2+1,$3,$7,$8,$9,$6}'
done |  sort -k1,1 -k4,5n > $output

echo "End"
date
echo ""