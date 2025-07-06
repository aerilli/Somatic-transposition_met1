#!/bin/bash

## 


source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/recombination_EVD

insertion_list="somatic-insertions.txt"

cd /ebio/abt6_projects/met1_somatic_transpositions/data/check_partial-mapping-per-TE_REVISIONS

GENOME=/ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.genome

> random-stretches.bed
lengths=`seq -s ' ' 500 500 10000`
for L in $lengths; do
 bedtools random -g $GENOME -l $L -n 10000 >> random-stretches.bed
done

> random-stretches_cx.bed
> random-stretches_ALL.bed
SP=`echo met1_0{1..9} met1_10 WT`
for s in $SP; do 
    BAM=/ebio/abt6_projects/met1_somatic_transpositions/data/2_met1s/2_aln/merged-runs/${s}.Chr.sorted.bam 
    bedtools coverage -a random-stretches.bed -b $BAM -f 0.9999 -counts | awk -v S="${s}" '{print $0 "\t" S "\t" "cx"}' >> random-stretches_cx.bed
    bedtools coverage -a random-stretches.bed -b $BAM -counts | awk -v S="${s}" '{print $0 "\t" S "\t" "rest"}' >> random-stretches_ALL.bed
done


cat random-stretches_cx.bed random-stretches_ALL.bed  | sortBed -i stdin > random-stretches_cx-ALL.bed && rm random-stretches_cx.bed random-stretches_ALL.bed

conda deactivate