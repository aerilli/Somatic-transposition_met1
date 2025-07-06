#!/bin/bash


#

# set -xv


cd /ebio/abt6_projects/met1_somatic_transpositions/data/RNA-Seq/pedestrian

> intermediate.libsizes-uniquely-mapped.tsv
> intermediate.libsizes-raw-count.tsv
> intermediate.ALL.intersect.counts.tsv
> ALL.weighted-counts.tsv

for i in  Tsu_P2_{1..3} Tsu_P1_{1..3} Tsu_WT_{1..3} Tsu_T4_P2_{1..3}; do

  echo ""
  echo Now doing $i

  if [ ! -f $i.ATCOPIA93-intersect.strict-bt2.sorted.bam ]; then 
  echo $i ; zcat reads/${i}*.fq.gz |
  bowtie2  --mp 13 --rdg 8,5 --rfg 8,5 --very-sensitive  --threads 32 \
  -x /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.scaffolds_contigs.v2  -U - |
  samtools sort -@24  -O BAM - | bedtools intersect -a stdin -b ATCOPIA93.Copia_LTR_retrotransposon.intact-EDTA.gff3 -wb   \
  > $i.ATCOPIA93-intersect.strict-bt2.sorted.bam 
  samtools index -@24 $i.ATCOPIA93-intersect.strict-bt2.sorted.bam
  fi


  # > libsizes.tsv
  ## uniquely mapped reads
  l=`zcat reads/${i}*.fq.gz |
  bowtie2  --mp 13 --rdg 8,5 --rfg 8,5 --very-sensitive  --threads 32 \
  -x /ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2 -U - |
  samtools sort -@32 |
  samtools view -@32 -f 4 |
  cut -f1 | sort | uniq | wc -l `

  echo $l >> intermediate.libsizes-uniquely-mapped.tsv

  ## All reads
  r=`zcat reads/${i}*.fq.gz | sort | wc -l `
  echo $r | tee >> intermediate.libsizes-raw-count.tsv

  # weighted count
  bedtools bamtobed -i $i.ATCOPIA93-intersect.strict-bt2.sorted.bam |
  awk '{OFS="\t"}{print $1,$2,$3,$5/60}' |
  bedtools intersect -a stdin -b ATCOPIA93.Copia_LTR_retrotransposon.intact-EDTA.gff3 -wao |
  awk -v I="${i}" -v L="${l}" -v R="${r}" '{OFS="\t"}{print $1,$2,$3,$4,$5":"$8"-"$9,L, R,I}' \
  sed 's/\(.*\)_/\1\t/' \
    >> intermediate.weighted-counts.tsv



  # > ALL.intersect.counts.tsv
  bedtools intersect -a ATCOPIA93.Copia_LTR_retrotransposon.intact-EDTA.gff3 -b $i.ATCOPIA93-intersect.strict-bt2.sorted.bam -c |
  awk -v I="${i}" -v L="${l}" -v R="${r}"  '{OFS="\t"}{print $0, $10/L, L, $10/R, R, I}' >> intermediate.ALL.intersect.counts.tsv




  echo ""
  echo Finished with $i

done

echo Now removing intermediates

# sorted by what?
cat intermediate.libsizes-uniquely-mapped.tsv| sort | uniq > libsizes-uniquely-mapped.tsv #&& rm intermediate.libsizes-uniquely-mapped.tsv
cat intermediate.libsizes-raw-count.tsv| sort | uniq > libsizes-raw-count.tsv #&& rm intermediate.libsizes-raw-count.tsv
cat intermediate.ALL.intersect.counts.tsv | sort | uniq | sort -k 16 > ALL.intersect.counts.tsv #&& rm intermediate.ALL.intersect.counts.tsv

for i in `awk '{print $8"-"$9}' intermediate.weighted-counts.tsv | sort | uniq `; do
 genotype=`echo ${i} | cut -f1 -d"-"`
 biorep=`echo ${i} | cut -f2 -d"-"`
 awk -v G="${genotype}" -v B="${biorep}" '{OFS="\t"}{print $5,"0",$5,1,1,G,B}' intermediate.weighted-counts.tsv | sort | uniq |
  sed 's/:/\t/' | sed 's/-/\t/' >>  intermediate.weighted-counts.tsv
done
cat intermediate.weighted-counts.tsv | sort | uniq | sort -k 6 > ALL.weighted-counts.tsv && rm intermediate.weighted-counts.tsv


cat ALL.intersect.counts.tsv  | sed 's/[A-Za-z0-9;=_\/:]*Name=//g' | sed 's/;/\t/' |
 sed 's/[A-Za-z0-9;=_\/:]*ltr_identity=//g' | sed 's/;[A-Za-z0-9;=_\/:]*//g' |
 awk '{OFS="\t"}{print $1":"$4,$9,$10, $11, $12, $13, $14, $15, $16}'  | sed 's/\(.*\)_/\1\t/' > ALL.intersect.counts.mod-4R.tsv
