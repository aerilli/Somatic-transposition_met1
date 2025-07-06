#!/bin/bash

source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/pbcoretools_2023

set -xv

joutname=$1
CORES=$2

WD=/ebio/scratch/amovilli/pbsim/ccs
input_file_FP=`ls -Art /ebio/scratch/amovilli/pbsim/method_errhmm/sd_merged.simulated-insertions.${joutname}.sorted.bam | tail -1`   # sd_merged.simulated-insertions.1.sorted.bam
input_file=sd_merged.simulated-insertions.${joutname}.sorted.bam 

mkdir -p $WD
cd $WD

# Generate circular consensus sequences (CCS) from subreads in n chunks, into a directory of your choice:
# source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/ccs
# date
# ccs --all -j $CORES --report-file $WD/${input_file%.sorted.bam}.ccs_report.txt $input_file_FP $WD/${input_file%.sorted.bam}.ccs.bam
# date
# conda deactivate

# pbindex ${input_file%.sorted.bam}.ccs.bam
# ### pbindex creates a index file that enables random-access to PacBio-specific data in BAM files.


# echo "Extracting HiFi-only reads"
# extracthifi ${input_file%.sorted.bam}.ccs.bam ${input_file%.sorted.bam}.hifi.bam
# ###  extract HiFi reads (>= Q20) from full CCS reads.bam output

# pbindex ${input_file%.sorted.bam}.hifi.bam

# bam2fastq -o ${input_file%.sorted.bam}.hifi ${input_file%.sorted.bam}.hifi.bam && rm sd_merged.simulated-insertions.ccs* && rm sd_merged.simulated-insertions.hifi.bam 





input_hifi=`realpath ${input_file%.sorted.bam}.hifi.fastq.gz`


# rm `ls | grep -v .hifi.fastq.gz`

# set +xv

conda deactivate


source activate /ebio/abt6_projects/met1_somatic_transpositions/conda/so_ALNtransposition

# set -xv

WD=/ebio/scratch/amovilli/pbsim/aln
REF=/ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.fa

mkdir -p $WD
cd $WD

# winnowmap mapping reads to WT genome
meryl count k=15 output $WD/merylDB $REF
meryl print greater-than distinct=0.9998 $WD/merylDB > $WD/merylDB/repetitive_k15.txt

nice  winnowmap -W repetitive_k15.txt -Y -L  -t $CORES -ax map-pb ${REF} ${input_hifi} | samtools view -h -b -  > $WD/Simulated_reads_insertions-aln2-WT-Chr_winnowmap.bam   # -F2308 -e '![SA]' are removed!!
samtools sort -@ $CORES  -O bam -o $WD/Simulated_reads_insertions-aln2-WT-Chr_winnowmap.sorted.bam $WD/Simulated_reads_insertions-aln2-WT-Chr_winnowmap.bam


samtools addreplacerg -w -r "@RG\tID:simulated\tLB:simulated\tPL:sequelII_hifi_v3\tSM:simulated" \
		-@ $CORES -o $WD/Simulated_reads_insertions-aln2-WT-Chr_winnowmap.sortedRG.bam \
			$WD/Simulated_reads_insertions-aln2-WT-Chr_winnowmap.sorted.bam -O BAM &&
			samtools index $WD/Simulated_reads_insertions-aln2-WT-Chr_winnowmap.sortedRG.bam &&
			rm $WD/Simulated_reads_insertions-aln2-WT-Chr_winnowmap.bam && rm $WD/Simulated_reads_insertions-aln2-WT-Chr_winnowmap.sorted.bam



# set +xv

conda deactivate