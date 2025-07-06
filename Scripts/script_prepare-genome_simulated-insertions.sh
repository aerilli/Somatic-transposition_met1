#!/bin/bash

source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/gffutils

set -xv

WD=/ebio/scratch/amovilli/pbsim
FREPEATS=/ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/3_annotation_REPEATS/Chr_scaffolds/Tsu-0/Tsu-0.Repeats_merged.Nslop50.nocentr.gff ## Apr 2025 Added filtering of repeats

cd $WD

mkdir -p sample

## new design for simulation! 24-08-16

EVDfa=`tail -1  /tmp/global2/amovilli/TE_resources/EVD_1.fa`

# generating fasta
> ./sample/genome.sim-insertions.fa
for ((i=1; i<=100; i++)) ; do 
bioawk -c fastx -v I="${i}" '{print ">" $name "_"I "\n" $seq }' /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.fa \
  >> ./sample/genome.sim-insertions.fa
done

# generating associated genome file
bioawk -c fastx '{print $name, length($seq)}' ./sample/genome.sim-insertions.fa > ./sample/genome.sim-insertions.genome

# generating 500 insertions, avg 1 per chr
> ./sample/random-ins.bed
for ((i=1; i<=500; i++)) ; do 
 echo insertion no $i 
 > ./sample/single-random-ins.bed
 while true ; do ## Apr 2025 Added loop
    bedtools random -g ./sample/genome.sim-insertions.genome -l 0 -n 1 | bedtools intersect -v -a stdin -b $FREPEATS > ./sample/single-random-ins.bed ## Apr 2025 Added filtering of repeats
    if [ -s ./sample/single-random-ins.bed ]; then
      break
    else
      sleep 1
    fi
 done

  while IFS= read -r line ; do
    chr=$(echo "${line}" | cut -f1)
    pos=$(echo "${line}" | cut -f2)

    echo "${line}" | cut -f1
    echo "${line}" | cut -f2

    cat ./sample/genome.sim-insertions.fa |
     seqkit mutate -i ${pos}:${EVDfa} -s $chr --quiet > ./sample/tmp.fa 
    #  exit 1 
    cat ./sample/tmp.fa > ./sample/genome.sim-insertions.fa && rm ./sample/tmp.fa

    # exit 1 
    bioawk -c fastx '{print $name, length($seq)}' ./sample/genome.sim-insertions.fa  > ./sample/genome.sim-insertions.genome

  done < ./sample/single-random-ins.bed

  cat ./sample/single-random-ins.bed >> ./sample/random-ins.bed

done

sort -k1,1 -k2,2n ./sample/random-ins.bed | cut -f 1-3 > ./sample/genome.sim-insertions.bed && rm ./sample/random-ins.bed

samtools index ./sample/genome.sim-insertions.fa

gzip -f ./sample/genome.sim-insertions.fa

set +xv 

conda deactivate