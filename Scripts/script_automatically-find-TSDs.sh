#!/bin/bash


REF=/ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.Chr_scaffolds.v2.fa
seq=/ebio/abt6_projects/met1_somatic_transpositions/data/2_met1s/1_lima_deepvariant/*.dc_q20.fastq.gz
reseq=/ebio/abt6_projects/met1_somatic_transpositions/data/2_met1s/reseq_1_lima_deep/*r.q20.fastq.gz


source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/gffutils

cat /ebio/abt6_projects/met1_somatic_transpositions/data/automatic-TSD-retrieval/SomaticInsertions_CURATED-list.no-comm.mod.24-06-14.CONFIRMED.tsv |
 cut -f12 |
 tr ',' '\n' \
  > list.reads-ID.txt

## LONG STEP, not to be repeated if I want only to generate the TSDs.txt file
# > reads.fa
# for i in $seq ; do 
# seqkit grep -j 32 -f list.reads-ID.txt  $i  | seqkit fq2fa  \
#  >> reads.fa
# done 
# for i in $reseq ; do 
# seqkit grep -j 32 -f list.reads-ID.txt  $i  | seqkit fq2fa  \
#  >> reads.fa
# done 

sed -i 's/ bc=.*//' reads.fa

# #set -xv
#### method 1 ####
echo -e "#TEfam\tTEcopy\tTSDlen\treadID\tChrinssite\tstartinssite\tendinssite" > TSDs.txt
while IFS= read -r line; do
    
    inssite=$(echo "${line}" | cut -f1,2,3)
    read=$(echo "${line}" | cut -f12)
    TEcopy=$(echo "${line}" | cut -f11)
    TEfam=$(echo "${line}" | cut -f10)
    
    for i in `echo $read | tr ',' '\t' `; do
        FirstSize=$(du -b TSDs.txt | cut -f 1)
        > inssite.fa
        > read.fa
        seqkit grep -p $i reads.fa > read.fa
        lenread=`grep -v "^>" read.fa | wc -c`
        
        # echo $read

        echo $inssite | tr ' ' '\t' |
        bedtools slop -i stdin -g ${REF%.fa*}.genome -b $lenread |
        bedtools getfasta -fi $REF -bed stdin \
        > inssite.fa

        # April 2025, after revision: blast to bed is a 1-based to 0-based conversion, therefore <awk '{OFS="\t"}{ if ($9<$10){print $2,$9,$10,$1} else if($9>$10){print $2,$10,$9,$1}}'> ----> <awk '{OFS="\t"}{ if ($9<$10){print $2,$9-1,$10,$1} else if($9>$10){print $2,$10-1,$9,$1}}'>
        > tmp.blastout.bed
        blastn -subject inssite.fa \
        -query read.fa \
        -outfmt 6 |
        awk '{if ($12>2000){print $0}}'| sort -k12,12 | tail -2 |
        awk '{OFS="\t"}{ if ($9<$10){print $2,$9-1,$10,$1} else if($9>$10){print $2,$10-1,$9,$1}}' |
        sortBed -i - |
        tee tmp.blastout.bed |
        bedtools merge -i stdin |
        bedtools coverage -a stdin -b tmp.blastout.bed -hist |
        awk -v I="${inssite}" -v R="${i}" -v C="${TEcopy}" -v T="${TEfam}" '{OFS="\t"}{if ($4==2 && $5<15){print T,C, $5, R, I} }' >> TSDs.txt
        

        if [ ! -s read.fa ] | [ ! -s tmp.blastout.bed ]; then
         echo -e "$TEfam\t$TEcopy\tnot-possible-to-retrieve\t$i\t$inssite"   >> TSDs.txt
        fi
        
        SecondSize=$(du -b TSDs.txt | cut -f 1)
        if [[ "$FirstSize" == "$SecondSize" ]]; then
         echo -e "$TEfam\t$TEcopy\tno-TSD-possibly-no-central-config\t$i\t$inssite" >> TSDs.txt 
        fi

        # echo "#" >> TSDs.txt
    done


done < /ebio/abt6_projects/met1_somatic_transpositions/data/automatic-TSD-retrieval/SomaticInsertions_CURATED-list.no-comm.mod.24-06-14.CONFIRMED.tsv




conda deactivate