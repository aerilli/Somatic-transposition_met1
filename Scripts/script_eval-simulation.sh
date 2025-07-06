#!/bin/bash

WD=$1    # /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/2024-09-20_cov2_sample

mkdir -p $WD/eval

cd $WD/eval

cp /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sim-insertions.bed ./sim-insertions.bed

source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/gffutils


GT=$WD/eval/sim-insertions.bed # groundtruth

> retrieved-simins.txt

for i in `seq -s ' ' 1 10 `; do # Change to how many iterations are done so far

    for j in `ls ../*.somatic-calls_CG.bed | sed -E 's@\.\.\/[0-9]+.@@ ; s/\.somatic.*//' | sort | uniq ` ; do

    > ${i}.${j}.somatic-calls_CG-SA.formatted.bed

    jmod=`echo $j | tr '-' '.'`

    cg=../${i}.${j}.somatic-calls_CG.bed
    sa=../${i}.${j}.somatic-calls_SA.bed

    cat $sa  | sed -e 's/\# DIFFICULT CASE TO AUTOMATE: S[0-9]*\/[0-9]*\/ccs [0-9]* [0-9]* Chr/Chr/' |tr ' ' '\t' | cut -f1-3 |
     awk -v JMOD="${jmod}"  -v I="${i}" '{OFS="\t"}{print $0,JMOD, I, "SA"}' >> ${i}.${j}.somatic-calls_CG-SA.formatted.bed
    tail +2 $cg | sed -e 's/\# DIFFICULT CASE TO AUTOMATE: S[0-9]*\/[0-9]*\/ccs [0-9]* [0-9]* Chr/Chr/' |tr ' ' '\t' | cut -f1-3 |
     awk -v JMOD="${jmod}"  -v I="${i}" '{OFS="\t"}{print $0,JMOD, I, "CG"}' >> ${i}.${j}.somatic-calls_CG-SA.formatted.bed
    
    sortBed -i ${i}.${j}.somatic-calls_CG-SA.formatted.bed > ${i}.${j}.somatic-calls_CG-SA.formatted.sorted.bed
    
    sed 's/_[0-9]*//' $GT |bedtools slop -i stdin -g /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.genome  -b 5 | 
     bedtools intersect -a stdin -b ${i}.${j}.somatic-calls_CG-SA.formatted.sorted.bed  | wc -l | awk -v JMOD="${jmod}"  -v I="${i}" '{OFS="\t"}{print $0, JMOD, I}' >> retrieved-simins.txt
    


    done
done

conda deactivate