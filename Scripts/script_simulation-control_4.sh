#!/bin/bash

## new design for simulation! 24-08-16

source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/gffutils


date
echo "Start!"
echo "Preparing genome..."

CORES=${1:-${NSLOTS:-32}}

# nice /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/script_prepare-genome_simulated-insertions.sh

# cp /ebio/scratch/amovilli/pbsim/sample/genome.sim-insertions.bed /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/sim-insertions.bed #&& rm /ebio/scratch/amovilli/pbsim/sample/genome.sim-insertions.bed


callsomatic() {
  i=$1
  joutname=$2

  echo "${i} ${j} SA calling..."
  

  /ebio/abt6_projects/met1_somatic_transpositions/code/scripts_re/script_re-calling_intact-TE_insertions_strict_2.sh \
    -s  ${i} \
    -b  /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/${i}.${joutname}.bam \
    -o  /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/somatic-insertions_call \
    -t  /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/5_TEanno/div_5_Col-CC_v2-2_EDTA/Tsu-0.scaffolds_contigs.v2.fa.mod.EDTA.TEanno.copy.merged-overlapping-elements.manual.24-05-30.gff \
    -f  /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/3_annotation_REPEATS/Chr_scaffolds/Tsu-0/Tsu-0.Repeats_merged.gff \
    -R  /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.fa \
     2>&1 | tee -a /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/logs/SA-ins_$(date '+%Y-%m-%d_%H')_${i}.${joutname}.log

  echo "${i} ${j} CG calling..."

  /ebio/abt6_projects/met1_somatic_transpositions/code/scripts_re/script_insertions_CIGAR.sh \
    -s  ${i} \
    -b  /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/${i}.${joutname}.bam \
    -o  /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/somatic-insertions_call \
    -t  /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/5_TEanno/div_5_Col-CC_v2-2_EDTA/Tsu-0.scaffolds_contigs.v2.fa.mod.EDTA.TEanno.copy.merged-overlapping-elements.manual.24-05-30.gff \
    -f  /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/3_annotation_REPEATS/Chr_scaffolds/Tsu-0/Tsu-0.Repeats_merged.gff \
    -R  /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.fa  \
    2>&1 | tee -a /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/logs/CG-ins_$(date '+%Y-%m-%d_%H')_${i}.${joutname}.log


  echo "${i} ${j} copying somatic calls..."
    cp  /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/somatic-insertions_call/${i}.${joutname}.TE-overlap.SA.TE-format.re-overlap.merged.final.bed.insertions.cov.bed /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/${i}.${joutname}.somatic-calls_SA.bed
    cp  /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/somatic-insertions_call/${i}.${joutname}.CG.insertions.final.bed /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/${i}.${joutname}.somatic-calls_CG.bed
    # rm /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/${i}.${joutname}.bam 


}





mkdir -p /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/

for i in `seq -s ' ' 4 10`; do #1 to 10

    j=1.0

    joutname=`echo $j | sed 's/\./-/'`

    date

    echo "Iteration number ${i}" | tee >(mail -s "Iteration number ${i} is starting!" "amovilli")

    echo "${i} ${j} hifi reads simulation..."

    /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/script_pbsim_3.sh ${j} $CORES

    date
    echo "${i} ${j} ccs and alignment..."

    /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/script_simulated-insertion_reads_ccs_aln_2.sh ${joutname} $CORES

    cp /ebio/scratch/amovilli/pbsim/aln/Simulated_reads_insertions-aln2-WT-Chr_winnowmap.sortedRG.bam /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/${i}.${joutname}.bam # && rm /ebio/scratch/amovilli/pbsim/aln/Simulated_reads_insertions-aln2-WT-Chr_winnowmap.sortedRG.bam*
    samtools index /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/${i}.${joutname}.bam || exit 1
    
    bam1=/ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/${i}.${joutname}.bam

    date
    echo "${i} ${j} somatic insertion calling..."

    # 1
    callsomatic $i $joutname


    
    # 0.1
    j=0.5
    joutname=`echo $j | sed 's/\./-/'`

    echo "${i} ${j} subsample reads..."
    date

    samtools view -@${CORES} -h -b --subsample $j  $bam1  |
     samtools sort -@${CORES} -O BAM - >  /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/${i}.${joutname}.bam
     samtools index /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/${i}.${joutname}.bam
    
    echo "${i} ${j} somatic insertion calling..."

    callsomatic $i $joutname


    # 0.01
    j=0.25
    joutname=`echo $j | sed 's/\./-/'`
    
    echo "${i} ${j} subsample reads..."
    date

    samtools view -@${CORES} -h --subsample $j  $bam1 -b  |
     samtools sort -@${CORES} -O BAM - >  /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/${i}.${joutname}.bam
    samtools index /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/${i}.${joutname}.bam

    echo "${i} ${j} somatic insertion calling..."
    
    callsomatic $i $joutname



    # 0.01
    j=0.05
    joutname=`echo $j | sed 's/\./-/'`
    date
    
    echo "${i} ${j} subsample reads..."

    samtools view -@${CORES} -h --subsample $j  $bam1 -b  |
     samtools sort -@${CORES} -O BAM - >  /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/${i}.${joutname}.bam
    samtools index /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/${i}.${joutname}.bam

    echo "${i} ${j} somatic insertion calling..."
    
    callsomatic $i $joutname
  

  # 0.001
    j=0.005
    joutname=`echo $j | sed 's/\./-/'`
    date
    
    echo "${i} ${j} subsample reads..."

    samtools view -@${CORES} -h --subsample $j  $bam1 -b  |
     samtools sort -@${CORES} -O BAM - >  /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/${i}.${joutname}.bam
    samtools index /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/${i}.${joutname}.bam

    echo "${i} ${j} somatic insertion calling..."
    
    callsomatic $i $joutname

# exit 1

    rm /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/sample/${i}.*.bam 


done


conda deactivate

date
echo "End!"