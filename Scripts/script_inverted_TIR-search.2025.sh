#!/bin/bash


cd /ebio/abt6_projects/met1_somatic_transpositions/data/inverted_TIR-search_REVISIONS

source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/gffutils

TEanno=/ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/5_TEanno/div_5_Col-CC_v2-2_EDTA/Tsu-0.scaffolds_contigs.v2.fa.mod.EDTA.TEanno.copy.merged-overlapping-elements.manual.24-06-13.gff
FREPEATS=/ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/3_annotation_REPEATS/Chr_scaffolds/Tsu-0/Tsu-0.Repeats_merged.gff


>close-TEs.same-family.gff

for i in `bedtools intersect -v  -a $TEanno -b $FREPEATS  | cut -f9 | sed 's/ID=.*Name=//;s/\;.*//' | sort | uniq` ; do
    
    echo $i
    grep $i $TEanno > unique_TE-family.gff
    > unique_TE-family.close.gff

    while IFS= read -r line; do
        #echo ""
        #echo $line
        echo "${line}" | bedtools closest -S -io -d -a stdin -b unique_TE-family.gff |
         awk '$NF <= 500 && $NF != -1' |
         sort -k1,1 -k4,4n |
         bedtools intersect -v  -a stdin -b $FREPEATS  |
         awk -v I="${i}" '{print $0 "\t" I}' >> unique_TE-family.close.gff
    done < unique_TE-family.gff

         bedtools complement -L -i  unique_TE-family.close.gff -g /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.genome | #bedtools merge -i stdin -d 500 |
         awk -v I="${i}"  '{if ($3-$2<500){print $0  "\t" I}}' >>close-TEs.same-family.gff

done


SP=`echo met1_0{1..9} met1_10 WT`
> close-TEs.same-family.2048-overlap.gff
for s in $SP ; do 

    echo $s

    BAM=/ebio/abt6_projects/met1_somatic_transpositions/data/2_met1s/2_aln/merged-runs/${s}.Chr.sorted.bam 
    samtools view -@64 -h -q 50 -f 2048 $BAM | cut -f1 | sort | uniq > $(basename ${BAM}  ".bam").2048.list
    samtools view -@64 -h -q 50 -N $(basename ${BAM}  ".bam").2048.list $BAM |bamToBed -i stdin > $(basename ${BAM}  ".bam").2048.bed

    bedtools intersect -v  -a $(basename ${BAM}  ".bam").2048.bed -b $FREPEATS |
    bedtools coverage -a close-TEs.same-family.gff    -b stdin -counts | awk -v S="${s}" '$5>0{print $0 "\t" S}' >> close-TEs.same-family.2048-overlap.gff

done

bedtools igv -i close-TEs.same-family.2048-overlap.gff > igv.batch

cp close-TEs.same-family.2048-overlap.gff close-TEs.same-family.2048-overlap.COMMENTS.gff 


realpath close-TEs.same-family.2048-overlap.COMMENTS.gff 