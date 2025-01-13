#!/bin/bash

# Virtual env activation:
# containing bioawk; bedtools v2.30.0; samtools v1.16.1
## source activate /ebio/abt6_projects/met1_somatic_transpositions/conda/so_somatic-events

mkdir ./tmp/SA_DEL
cd ./tmp/SA_DEL

samtools merge -@32  -o ALL_met1.bam `find  ./Data/ALN -name "met1_*.Chr.sorted.bam"`

samtools view -@32 -h ALL_met1.bam -f 2048 | cut -f1 > SA.list
samtools view -@32 -hb ALL_met1.bam -N SA.list  | samtools sort -@32 -O BAM - > ALL_met1.SA.bam

samtools index ALL_met1.SA.bam

samtools merge -o ALL_met1.DEL.bam `find ./Data/Somatic_events/excision_CIGAR -maxdepth 1 -name "met1_*.Chr.sorted.excision.500.bam"`


bedtools intersect -a ALL_met1.SA.bam -b ALL_met1.DEL.bam -bed |
 awk '{if ($5>50){print $0}}'| sortBed |
  bedtools merge -i stdin | tee ALL_met1.SA-DEL.bed |
  bedtools intersect -v -a stdin -b ./Data/Tsu-0.repeats.merged.gff \
    > ALL_met1.SA-DEL.no-rep.bed


## conda deactivate