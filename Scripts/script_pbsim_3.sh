#!/bin/bash

source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/pbsim3

set -xv

subs=$1 # for subsampling simulation reads.
CORES=${2:-128}

subsoutname=`echo $subs | sed 's/\./-/'`
WD=/ebio/scratch/amovilli/pbsim

cd $WD

mkdir -p sample

mkdir -p $WD/method_errhmm

#gzip  -d ./sample/genome.sim-insertions.fa.gz

cd ${WD}/method_errhmm

date
pbsim --strategy wgs \
      --method errhmm \
      --length-mean 15775.0 \
      --length-sd 7500 \
      --difference-ratio 22:45:33 \
      --length-min 5000 \
      --errhmm /ebio/scratch/amovilli/pbsim/pbsim3/data/ERRHMM-SEQUEL.model \
      --depth 4 \
      --genome ../sample/genome.sim-insertions.fa \
      --accuracy-mean 0.9 \
      --pass-num 7
date
set +xv

conda deactivate

# gzip --fast -f ../sample/genome.sim-insertions.fa

source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/pbcoretools_2023


set -xv

echo sorting $subs
date

mkdir -p /ebio/scratch/amovilli/pbsim/ccs/logs

for i in sd_0*.sam ; do 
  # samtools view -@ ${CORES} -hbS $i | samtools sort -@ ${CORES} -O BAM - > ${i%.sam}.${subsoutname}.sorted.bam ;
  # samtools index ${i%.sam}.${subsoutname}.sorted.bam ;
  # ccs --all -j $CORES --report-file  /ebio/scratch/amovilli/pbsim/ccs/${i%.sam}.${subsoutname}.ccs_report.txt ${i%.sam}.${subsoutname}.sorted.bam  /ebio/scratch/amovilli/pbsim/ccs/${i%.sam}.${subsoutname}.sorted.ccs.bam 
  # pbindex /ebio/scratch/amovilli/pbsim/ccs/${i%.sam}.${subsoutname}.sorted.ccs.bam 
  # extracthifi /ebio/scratch/amovilli/pbsim/ccs/${i%.sam}.${subsoutname}.sorted.ccs.bam  /ebio/scratch/amovilli/pbsim/ccs/${i%.sam}.${subsoutname}.sorted.hifi.bam 
  # pbindex /ebio/scratch/amovilli/pbsim/ccs/${i%.sam}.${subsoutname}.sorted.hifi.bam 	
	
	qsub -pe parallel 4 -N ${i%.*} -l h_vmem="10G" \
   -o /ebio/scratch/amovilli/pbsim/ccs/logs/$(date '+%Y-%m-%d-%H').ccs_${i%.*}.log -j y \
   /ebio/abt6_projects/Ath_HiFi_met1/code/submit_eNO.bash --submit_command \
   "$(echo /ebio/abt6_projects/met1_somatic_transpositions/data/insertions-simulation/script_ccs-submit.sh "${i}" "${subsoutname}" 4 )"

	
done

while qstat | awk 'NR > 2 && $5 != "Eqw"' | grep -q .; do
    echo "Jobs still running... waiting"
    sleep 60  # check every 30 seconds
done

echo "All jobs finished!" | mail -s "ccs succesful"  amovilli
  



date

pbmerge /ebio/scratch/amovilli/pbsim/ccs/*.${subsoutname}.sorted.hifi.bam > /ebio/scratch/amovilli/pbsim/ccs/sd_merged.simulated-insertions.${subsoutname}.hifi.bam
pbindex /ebio/scratch/amovilli/pbsim/ccs/sd_merged.simulated-insertions.${subsoutname}.hifi.bam
bam2fastq -o /ebio/scratch/amovilli/pbsim/ccs/sd_merged.simulated-insertions.${subsoutname}.hifi /ebio/scratch/amovilli/pbsim/ccs/sd_merged.simulated-insertions.${subsoutname}.hifi.bam  && rm sd_0* 
date
# set -xv
# echo "merging $subs into sd_merged.simulated-insertions.${subsoutname}.sorted.bam"

# # samtools merge -@ ${CORES} -f -O BAM sd_merged.simulated-insertions.${subsoutname}.bam  sd_0*.${subsoutname}.sorted.bam 
# #  samtools sort -@ $CORES  -O BAM -o sd_merged.simulated-insertions.${subsoutname}.sorted.bam   sd_merged.simulated-insertions.${subsoutname}.bam   && rm sd_merged.simulated-insertions.${subsoutname}.bam  && rm sd_0*.${subsoutname}.sorted.bam && echo "merging successful"
# # samtools index sd_merged.simulated-insertions.${subsoutname}.sorted.bam && rm sd_0*  && echo "indexing and clean-up successful"

# date
# ls sd_0*.${subsoutname}.sorted.bam  >files.bamlist
# bamtools merge -list files.bamlist -out sd_merged.simulated-insertions.${subsoutname}.bam
# date
# samtools sort -@ $CORES  -O BAM -o sd_merged.simulated-insertions.${subsoutname}.sorted.bam sd_merged.simulated-insertions.${subsoutname}.bam #&& rm sd_merged.simulated-insertions.${subsoutname}.bam  && rm sd_0*.${subsoutname}.sorted.bam && echo "merging successful"
# date

# samtools view -H sd_merged.simulated-insertions.${subsoutname}.sorted.bam |
# awk '
# /^@RG/ {
#     orig = $0
#     for (i = 1; i <= 500; i++) {
#         line = orig
#         sub(/PU:S1/, "PU:S" i, line)
#         print line
#     }
#     next
# }
# { print }
# ' > new_header.sam

# samtools reheader new_header.sam sd_merged.simulated-insertions.${subsoutname}.sorted.bam > sd_merged.simulated-insertions.${subsoutname}.sorted.tmp.bam  

# mv sd_merged.simulated-insertions.${subsoutname}.sorted.tmp.bam  sd_merged.simulated-insertions.${subsoutname}.sorted.bam  

# samtools index sd_merged.simulated-insertions.${subsoutname}.sorted.bam #&& rm sd_0*  && echo "indexing and clean-up successful"
date

set +xv






conda deactivate