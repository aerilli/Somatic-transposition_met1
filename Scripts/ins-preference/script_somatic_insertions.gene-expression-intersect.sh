#!/bin/bash

# set -xv


WD="/ebio/abt6_projects/met1_somatic_transpositions/data/3_somatic-events/r_mutations/insertion-bias/Expression_bias"

cd $WD

# gene list in organellar and other repeats
geneanno=/ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/6_geneAnno/gene-anno/Tsu-0/2407_liftoff/all-iso/Tsu-0.TAIR10.all-iso.liftoff.gff3
FREPEATS=/ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/3_annotation_REPEATS/Chr_scaffolds/Tsu-0/Tsu-0.Repeats_merged.gff

bedtools intersect -v -a $geneanno -b $FREPEATS |
 awk '{if ($3=="gene"){print $9}}' |
 sed 's/ID=//;s/;.*//' |
 grep -v ATM | grep -v ATC \
  > ./list-genes_not-overlapping-repeats.txt

grep -f ./list-genes_not-overlapping-repeats.txt ./rsem_25_07_2024_Tsu-0_run186_TPM.txt > ./rsem_25_07_2024_Tsu-0_run186_TPM.filtered.txt 



grep CONFIRMED ../SomaticInsertions_CURATED-list.no-comm.mod.tsv |
 cut -f1-3,10,11,13 |
 bedtools intersect -a stdin -b $geneanno -loj |
 awk '{OFS="\t"}{ if ($9=="gene") {print $1,$2,$3,$4,$5,$6, $15}}' |
 sed 's/ID=//;s/;.*//' \
  > somatic_insertions.gene-expression-intersect.bed



bedtools random -l 0 -n 10000 -g /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.genome -seed 888 |
 bedtools intersect -a stdin -b $geneanno -loj |
 awk -v R="Random" '{OFS="\t"}{ if ($9=="gene") {print $1,$2,$3,R,R,R, $15}}' |
 sed 's/ID=//;s/;.*//' \
  > rand_insertions.gene-expression-intersect.bed