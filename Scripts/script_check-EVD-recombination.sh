#!/bin/bash

source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/recombination_EVD
# mamba install -c bioconda minimap seqkit bcftools bedtools

set -xv

cd /ebio/abt6_projects/met1_somatic_transpositions/data/EVD_recombination

insertion_list="somatic-insertions.txt"


minimap2 -ax asm5 -t 4 --eqx EVD_Col-0.fa  EVD_1.fa > 1-to-col.sam
minimap2 -ax asm5 -t 4 --eqx EVD_Col-0.fa  EVD_2.fa > 2-to-col.sam

# cat EVD_1.fa EVD_2.fa | minimap2 -ax asm5 -t 4 --eqx EVD_Col-0.fa  -  > 1-2-to-col.sam
# bcftools mpileup -Ob -o 1-2-to-col.bcf -f EVD_Col-0.fa 1-to-col.sam 2-to-col.sam
# bcftools call -vmO v -o 1-2-to-col.vcf 1-2-to-col.bcf


bcftools mpileup -Ob -o 1-to-col.bcf -f EVD_Col-0.fa 1-to-col.sam 
bcftools call -vmO v -o 1-to-col.vcf 1-to-col.bcf

bcftools mpileup -Ob -o 2-to-col.bcf -f EVD_Col-0.fa 2-to-col.sam
bcftools call -vmO v -o 2-to-col.vcf 2-to-col.bcf

bedtools intersect -header -v  -a 1-to-col.vcf -b 2-to-col.vcf  > 1-to-col.no-2.vcf 
bgzip -f -k 1-to-col.no-2.vcf 
bcftools index 1-to-col.no-2.vcf.gz

bedtools intersect -header -v  -b 1-to-col.vcf -a 2-to-col.vcf  > 2-to-col.no-1.vcf 
bgzip -f -k 2-to-col.no-1.vcf 
bcftools index 2-to-col.no-1.vcf.gz

counter=0
grep -v "#" 1-to-col.no-2.vcf  | while read p; do
  ((counter++))
  grep "#" 1-to-col.no-2.vcf > ${counter}.SNP_EVD_1.vcf
  echo "$p" >> ${counter}.SNP_EVD_1.vcf
  bgzip -f -k ${counter}.SNP_EVD_1.vcf
  bcftools index ${counter}.SNP_EVD_1.vcf.gz
done 

counter=0
grep -v "#" 2-to-col.no-1.vcf  | while read p; do
  ((counter++))
  grep "#" 2-to-col.no-1.vcf > ${counter}.SNP_EVD_2.vcf
  echo "$p" >> ${counter}.SNP_EVD_2.vcf
  bgzip -f -k ${counter}.SNP_EVD_2.vcf
  bcftools index ${counter}.SNP_EVD_2.vcf.gz
done 



echo -e "#readid\tind\tinsertionas\t$(ls *.SNP_EVD_1.vcf | sed 's/.SNP_EVD_1.vcf/_1/'| sort -n | tr '\n' '\t')$(ls *.SNP_EVD_2.vcf | sed 's/.SNP_EVD_2.vcf/_2/'| sort -n |tr '\n' '\t')" > out.txt

for copy in 1 2 ; do
if [[ $copy == 1 ]]; then
grep "Evade" $insertion_list | grep "Chr1;11941" > reads_EVD_${copy}.txt
elif [[ $copy == 2 ]]; then
grep "Evade" $insertion_list | grep "Chr5;21419" > reads_EVD_${copy}.txt
fi

### how to treat commas
while IFS= read -r line; do
    insertionas="EVD_${copy}"
	  readids=$(echo "${line}" | cut -f8)
    ind=$(echo "${line}" | cut -f9 | sed 's/met1_0//;s/met1_//')


    for readid in `echo $readids | sed 's/,/\t/g'`; do
    
    echo -en "$readid\t$ind\t$insertionas\t" >> out.txt
    
    readid_prefix=`echo $readid | sed 's/\//\.\./g'`
    zgrep -q "${readid}" /ebio/abt6_projects/small_projects/amovilli/Tsu-0v2/1_lima/*_${ind}.dc_q20.fastq.gz &&
      seqkit grep -p "${readid}" /ebio/abt6_projects/small_projects/amovilli/Tsu-0v2/1_lima/*_${ind}.dc_q20.fastq.gz |
      minimap2 -ax asm5 -t 4 --eqx EVD_Col-0.fa  - > ${readid_prefix}-to-col.sam ||
        seqkit grep -p "${readid}" /ebio/abt6_projects/small_projects/amovilli/met1_Tsu-0_re-run/*_${ind}r.dc_q20.fastq.gz |
        minimap2 -ax asm5 -t 4 --eqx EVD_Col-0.fa  - > ${readid_prefix}-to-col.sam
      
      bcftools mpileup -Ob -o  ${readid_prefix}-to-col.bcf -f EVD_Col-0.fa  ${readid_prefix}-to-col.sam  ${readid_prefix}-to-col.sam
      bcftools call -vmO v -o  ${readid_prefix}-to-col.vcf  ${readid_prefix}-to-col.bcf
      bgzip -f -k ${readid_prefix}-to-col.vcf 
      bcftools index ${readid_prefix}-to-col.vcf.gz

      for i in `ls *.SNP_EVD_1.vcf | sed 's/.SNP_EVD_1.vcf//' | sort `; do
          
          io=`bcftools isec -n~11  ${readid_prefix}-to-col.vcf.gz  ${i}.SNP_EVD_1.vcf.gz | wc -l`
          echo -en "${io}\t" >> out.txt

      done

      for i in `ls *.SNP_EVD_2.vcf | sed 's/.SNP_EVD_2.vcf//' | sort `; do
          
          io=`bcftools isec -n~11 ${readid_prefix}-to-col.vcf.gz  ${i}.SNP_EVD_2.vcf.gz | wc -l`
          echo -en "${io}\t" >> out.txt

      done
      echo "" >> out.txt
    done

    # match1=`bedtools intersect -a $readid-to-col.vcf -b 1-to-col.no-2.vcf | wc -l`
    # match2=`bedtools intersect -a $readid-to-col.vcf -b 2-to-col.no-1.vcf | wc -l`





done < reads_EVD_${copy}.txt
done


set +xv


conda deactivate