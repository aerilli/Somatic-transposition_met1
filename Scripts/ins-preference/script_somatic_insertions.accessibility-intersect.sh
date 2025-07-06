#!/bin/bash

set -xv


WD="/ebio/abt6_projects/met1_somatic_transpositions/data/3_somatic-events/r_mutations/insertion-bias/Accessibility_bias"

cd $WD




for i in Tsu_P2_*_ref_Tsu0_P2gDNA_q005_SchmitzParam_treat_pileup.bdg ; do

    j=${i%_ref_Tsu0_P2gDNA_q005_SchmitzParam_treat_pileup.bdg}


    # if TE fell in btw annotations groupby > mean
    grep CONFIRMED  ../SomaticInsertions_CURATED-list.no-comm.mod.tsv |  cut -f1-3,10,11,13 |
     bedtools slop -i stdin -b 25 -g /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.genome |
     bedtools intersect -a stdin -b ${i}  -loj |
     cut -f 1-6,10 |
     bedtools groupby -i stdin -g 1-6 -c 7 -o mean |
     awk -v J="${j}" '{print $0 "\t" J}' \
      > ${j}.somatic_insertions.bslop25bp.accessibility-intersect.bed
     
 bedtools random -l 25 -n 10000 -g /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.genome -seed 888 |
     bedtools intersect -a stdin -b ${i}  -loj |
     awk -v R="Random" '{OFS="\t"}{print  $1,$2,$3,R,R,R,$10}' |
     bedtools groupby -i stdin -g 1-6 -c 7 -o mean |
     awk -v J="${j}" '{print $0 "\t" J}' \
      > ${j}.10kRandom_insertions.bslop25bp.accessibility-intersect.bed

done


