#!/bin/bash

# amovilli 2022-09-22

#### ONLY BATCH 1 ####


# USAGE # /ebio/abt6_projects9/Ath_HiFi_met1/code/mapping-scripts/ATAC-Seq/script_ATAC-Seq_mapping_merging-fq.sh <WORKDIR> <REFERENCE> <CORES>
# EXAMPLE # /ebio/abt6_projects9/Ath_HiFi_met1/code/mapping-scripts/ATAC-Seq/script_ATAC-Seq_mapping_merging-fq.sh  /ebio/abt6_projects9/Ath_HiFi_met1/data/mapping/BS-Seq/ Est-1 10



# Debugging
# set -xv

# Parse parameteres
WORKDIR=$1
ACC=$2
CORES=${3:-${NSLOTS:-6}}


# Declare other variables
output=$WORKDIR/$ACC/batch1
TSACC=`echo $ACC | sed 's/-//;s/[0-9]//g'` #remove - and all the numbers
if `echo $TSACC | grep -q "MAR"`
then
    TSACC="Mar"
fi

# Create necessary directories
mkdir -p /ebio/abt6_projects/met1_somatic_transpositions/data/3_somatic-events/r_mutations/insertion-bias/Accessibility_bias/ATAC-Seq/$ACC/batch1

mkdir -p $output/dataset/ && mkdir -p $output/genome/


echo ""
echo "fetching ATAC-Seq data"

cd $output/dataset


# # Activate virtual env
source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/RNA 


SCAFFOLD="/ebio/abt6_projects/met1_somatic_transpositions/data/1_genome/2_hifiasm_ragtag/Tsu-0/Tsu-0.scaffolds_contigs.v2.fa"
if test -f "$SCAFFOLD"; then
    ln -sf $SCAFFOLD $output/genome/$ACC.genome.fasta

fi


cd $output/genome

if [ `ls -1 ${ACC}*.bt* 2>/dev/null | wc -l ` -gt 0 ]; then
    echo "bt index is already present..."
else
    echo "Generating bt2 index for $ACC hifi genome..."
    bowtie2-build $output/genome/$ACC.genome.fasta $ACC.genome
fi



cd $output/dataset

while IFS= read -r line; do 

    if [[ $line =~ $TSACC ]] 
    then
    name=$(echo "${line}" | cut -f1)
    codname=$(echo "${line}" | cut -f2)

    R1=`printf "./""$codname"".R1.fastq.gz "`
    R2=`printf "./""$codname"".R2.fastq.gz "`
    R=`printf "./""$name.$codname"".R"`

    for n in {1..8}; do
        ATAC_R1=`printf "/ebio/abt6_projects/Bisulfite_RNA_ATAC_met1/data/ENAupload/ena_upload_2022-06-02_ATACseqBatch1/""$codname"".R1.RUN0190.L""$n"".fastq.gz "`
        ATAC_R2=`printf "/ebio/abt6_projects/Bisulfite_RNA_ATAC_met1/data/ENAupload/ena_upload_2022-06-02_ATACseqBatch1/""$codname"".R2.RUN0190.L""$n"".fastq.gz "`
        
        echo "Concatenating ${ATAC_R1} to ${R1}"
        if [ -f "$R1" ]; then echo "${R1} already present"; rm $R1; fi
        cat $ATAC_R1 >> $R1

        echo "Concatenating ${ATAC_R2} to ${R2}"
        if [ -f "$R2" ]; then echo "${R2} already present"; rm $R2; fi
        cat $ATAC_R2 >> $R2

    done



    echo "Trimming..."
    skewer -e -l 30 -L 100 -n -o $R -z -t $CORES $R1 $R2 && rm $R1 $R2



    echo "Aligning reads to $ACC hifi genome..."
    bowtie2 -x $output/genome/$ACC.genome -p $CORES -1 $R-trimmed-pair1.fastq.gz -2 $R-trimmed-pair2.fastq.gz  -S ${R}_mapped.sam
    samtools view -@ $CORES -q 2 -b ${R}_mapped.sam > ${R}_mapped.bam
    samtools sort -@ $CORES -l 9 ${R}_mapped.bam > ${R}_mapped.sorted.bam
    samtools view -@ $CORES ${R}_mapped.sorted.bam | cut -f2 | sort | uniq -c > ${R}_mapped.sorted.tagresults.txt

    samtools view -@ $CORES -b -q 30 -F 4 ${R}_mapped.sorted.bam > ${R}_mapped.sorted.filtered.bam 
    samtools view -@ $CORES -b -q 30 -F 268 ${R}_mapped.sorted.bam > ${R}_mapped.sorted.filtered.bam

    source activate /ebio/abt6_projects/Ath_HiFi_met1/conda/envs/picard
    picard MarkDuplicates I=${R}_mapped.sorted.filtered.bam \
    O=${R}_mapped.sorted.filtered.DEDUP.bam \
    M=${R}_metrics.txt \
    REMOVE_DUPLICATES=true

    conda deactivate

    samtools index ${R}_mapped.sorted.DEDUP.bam

    rm $R-trimmed-pair1.fastq.gz $R-trimmed-pair2.fastq.gz ${R}_mapped.sam ${R}_mapped.bam


    
    echo ""
    echo "=======================================================$name done======================================================="
    echo ""
   
    fi
   
 

done < /ebio/abt6_projects/met1_somatic_transpositions/data/3_somatic-events/r_mutations/insertion-bias/Accessibility_bias/ATAC-Seq/ATAC_batch1.Tsu_WT.txt



# # Deactivate virtual env
conda deactivate




echo ""
echo "=======================================================  DONE  ======================================================="
echo ""


 
# End
date
echo "End script"
echo ""
