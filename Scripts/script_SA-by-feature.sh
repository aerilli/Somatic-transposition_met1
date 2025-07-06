#!/bin/bash

## 


source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/recombination_EVD


set -xv
cd /ebio/abt6_projects/met1_somatic_transpositions/data/SA-by-feature_REVISIONS

awk '$3=="gene" {print $0}' "/ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/6_geneAnno/gene-anno/Tsu-0/2407_liftoff/single-iso/Tsu-0.TAIR10.single-iso.liftoff.gff3" > ./genes.gff
GENES="genes.gff"
grep rDNA /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/3_annotation_REPEATS/Chr_scaffolds/Tsu-0/Tsu-0.Repeats_merged.gff > ./rDNA.gff
rDNA=rDNA.gff
grep centromere /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/3_annotation_REPEATS/Chr_scaffolds/Tsu-0/Tsu-0.Repeats_merged.gff > ./centromere.gff
CEN=centromere.gff
grep telomere /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/3_annotation_REPEATS/Chr_scaffolds/Tsu-0/Tsu-0.Repeats_merged.gff > ./telomere.gff
TEL=telomere.gff
grep -e chloroplast -e mitochondria /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/3_annotation_REPEATS/Chr_scaffolds/Tsu-0/Tsu-0.Repeats_merged.gff > ./organellar.gff
ORG=organellar.gff
ln -s /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/5_TEanno/div_5_Col-CC_v2-2_EDTA/Tsu-0.scaffolds_contigs.v2.fa.mod.EDTA.TEanno.copy.merged-overlapping-elements.manual.24-06-13.gff ./TEs.gff
TES=TEs.gff

awk '{OFS="\t"}{print $1, 0, $2}' /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.genome > genome.bed
cat $TES /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/6_geneAnno/gene-anno/Tsu-0/2407_liftoff/single-iso/Tsu-0.TAIR10.single-iso.liftoff.gff3 \
 /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/3_annotation_REPEATS/Chr_scaffolds/Tsu-0/Tsu-0.Repeats_merged.gff |
 sortBed -i stdin | bedtools subtract  -a genome.bed -b stdin > ./intergenic.bed
INTERG=intergenic.bed

SP=`echo met1_0{1..9} met1_10 WT`
> ALL.intersect.2048.cov.bed
for s in $SP ; do 

    BAM=/ebio/abt6_projects/met1_somatic_transpositions/data/2_met1s/2_aln/merged-runs/${s}.Chr.sorted.bam 
    samtools view -@64 -h -q 50 -f 2048 $BAM | cut -f1 | sort | uniq > $(basename ${BAM}  ".bam").2048.list

    samtools view -@64 -h -q 50 -N $(basename ${BAM}  ".bam").2048.list $BAM |bamToBed -i stdin > $(basename ${BAM}  ".bam").2048.bed
    samtools view -@64 -h -q 50 $BAM |bamToBed -i stdin  > $(basename ${BAM}  ".bam").ALL.bed
    

    FEATS="$GENES $rDNA $CEN $TEL $ORG $TES $INTERG"
    for FEAT in $FEATS ; do

        bedtools merge -i $FEAT | bedtools coverage -a stdin -b $(basename ${BAM}  ".bam").2048.bed -counts | awk -v S="${s}" -v F="${FEAT%.*}" '{print $0 "\t" S "\t" F }' >> ALL.intersect.2048.cov.bed 
        bedtools merge -i $FEAT | bedtools coverage -a stdin -b $(basename ${BAM}  ".bam").ALL.bed -counts | awk -v S="${s}"  '{print $0 "\t" S "\t" "all" }'  >> ALL.intersect.2048.cov.bed 


    done
        rm  $(basename ${BAM}  ".bam").2048.bed $(basename ${BAM}  ".bam").ALL.bed
done 

set +xv

conda deactivate