#!/bin/bash

WORKDIR=$1 #/ebio/abt6_projects/met1_somatic_transpositions/data/BS-Seq

REF="/ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.fa"

mkdir -p $WORKDIR/dataset
mkdir -p $WORKDIR/genome

source activate /ebio/abt6_projects9/Ath_HiFi_met1/conda/envs/Bismark

cd $WORKDIR

ln -sf $REF `realpath $WORKDIR/genome/. `

bismark_genome_preparation --path_to_aligner /ebio/abt6_projects9/abt6_software/bin/bowtie2-2.2.3/ --verbose $WORKDIR/genome


PE=`find /ebio/abt6_projects/Bisulfite_RNA_ATAC_met1/data/ENAupload/ena_upload_2022-06-02_BSseqBatch1 -name "*Tsu0_P2*" -printf "%p "`
outname="Tsu0_P2"

cd dataset
skewer -q 20 -l 30 -n -z -o $WORKDIR/dataset/${outname} -t 32  ${PE}
bismark $WORKDIR/genome -1 $WORKDIR/dataset/${outname}-trimmed-pair1.fastq.gz -2 $WORKDIR/dataset/${outname}-trimmed-pair2.fastq.gz  && rm $WORKDIR/dataset/${outname}-trimmed-pair*.fastq.gz

deduplicate_bismark --bam -p $WORKDIR/dataset/${outname}-trimmed-pair1_bismark_bt2_pe.bam && rm $WORKDIR/dataset/${outname}-trimmed-pair1_bismark_bt2_pe.bam

cd ..
bismark_methylation_extractor -p  --cytosine_report --CX_context --genome_folder $WORKDIR/genome  $WORKDIR/dataset/${outname}-trimmed-pair1_bismark_bt2_pe.deduplicated.bam &&  rm $WORKDIR/${outname}-trimmed_bismark_bt2.deduplicated.bam

bismark2bedGraph -o Tsu0_P2.CpG CpG_OT_Tsu0_P2-trimmed-pair1_bismark_bt2_pe.deduplicated.txt CpG_OB_Tsu0_P2-trimmed-pair1_bismark_bt2_pe.deduplicated.txt && rm CpG_OT_Tsu0_P2-trimmed-pair1_bismark_bt2_pe.deduplicated.txt CpG_OB_Tsu0_P2-trimmed-pair1_bismark_bt2_pe.deduplicated.txt

bismark2bedGraph --CX_context -o Tsu0_P2.CHG CHG_OT_Tsu0_P2-trimmed-pair1_bismark_bt2_pe.deduplicated.txt CHG_OB_Tsu0_P2-trimmed-pair1_bismark_bt2_pe.deduplicated.txt && rm CHG_OT_Tsu0_P2-trimmed-pair1_bismark_bt2_pe.deduplicated.txt CHG_OB_Tsu0_P2-trimmed-pair1_bismark_bt2_pe.deduplicated.txt

bismark2bedGraph --CX_context -o Tsu0_P2.CHH CHH_OT_Tsu0_P2-trimmed-pair1_bismark_bt2_pe.deduplicated.txt CHH_OB_Tsu0_P2-trimmed-pair1_bismark_bt2_pe.deduplicated.txt && rm CHH_OT_Tsu0_P2-trimmed-pair1_bismark_bt2_pe.deduplicated.txt CHH_OB_Tsu0_P2-trimmed-pair1_bismark_bt2_pe.deduplicated.txt

conda deactivate