#!/bin/bash

# amovilli

## for detecting excised DNA TEs


set -xuv
# set -o pipefail


## EXAMPLE ###
## for SP in `cat ./list.SINGLE-PLANTS.txt` ; do
### ./Scripts/script_calling_somatic-excisions.CIGAR.sh\
###  -s  ${SP} \
###  -b  ./Data/ALN/${SP}-aln2-WT-Chr_winnowmap.sortedRG.bam \
###  -o  ./Data/Somatic_events/excisions_CIGAR \
###  -t  ./Data/Tsu-0.TE-anno.manual.gff \
###  -f  ./Data/Tsu-0.repeats.merged.gff \
###  -R  ./Reads/Tsu-0.reference.fa \
###   2>&1 | tee -a ./Data/logs/$(date '+%Y-%m-%d_%H-%M')_calling_TE_insertions_${SP}.log
## done


# setting variables
while getopts "s:b:o:t:f:R:h" opt; do
 case $opt in
s)
  singleplant=$OPTARG >&2
  ;;
b)
  BAMIN=$OPTARG >&2
  ;;
o)
  outdir=$OPTARG >&2
  ;;
t)
  TEanno=$OPTARG >&2
  ;;
f)
  FREPEATS=$OPTARG >&2
  ;;
R)
  REF=$OPTARG >&2
  ;;

h)
  echo "Help"
	exit 1
 esac
done


# Creating dirs and setting other variables 
mkdir -p $outdir
cd $outdir

BAMOUT=`basename ${BAMIN%.*}`
FREPEATSNS=${FREPEATS%.gff*}.Nslop50.nocentr.gff

# Virtual env activation:
# containing bioawk; bedtools v2.30.0; samtools v1.16.1
## source activate /ebio/abt6_projects/met1_somatic_transpositions/conda/so_somatic-events



#### #### ( 1 ) #### ####

# Increase N-stretch annotation of 50 bp in FREPEATS, beacuse variants are often found in proximity to GAPS. 
## Centromeres are retained because ATHILAs seem to mobilize (-> to be confirmed)
head -2 $FREPEATS \
  > $FREPEATSNS

bioawk -c gff '!/centromere/{if ( $feature == "N_stretch"){print $seqname, $source, $feature, $start-50, $end+50, $score, $filter, $strand, $group}  else {print $0}}'   $FREPEATS |
  grep Chr >> $FREPEATSNS






#### #### ( 2 ) #### ####
# from BAM input, extracting only reads that have a deletion ertion >500bp 
samtools view -@ 32 -h $BAMIN | awk 'BEGIN{OFS="\t"} { if ($6 ~ /[5-9][0-9][0-9]+D |[1-9][0-9][0-9][0-9]+D/) {print $0} else if ($1 ~ /^@/){print $0} }' |
 samtools view -@ 32 -q 50 -bS > $BAMOUT.excision.500.bam

samtools index $BAMOUT.excision.500.bam









#### #### ( 3 ) #### ####
# extracting read IDs with DEL  that overlap with TEanno
## at the same time I exclude those TE annotations that overlap with the adjusted (increased NStretches) FREPEATS
cat $TEanno | # Have to remove ChrM and ptg entries first
bedtools intersect -v -a - -b $FREPEATSNS | # removing TEs that overlap with FRepeats (rDNA, slop Nstretches, organellar, telomeres) from TE intact; centromeres included though
bedtools intersect -a $BAMOUT.excision.500.bam -b - -bed -wa -wb | cut -f4 | sort |uniq > read-IDs.reads.withExc.TE-overlap.list




############## DEL ##############

#### #### ( 4 ) #### ####
# filtering reads that have DEL  and overlap with intact TEs using the list above
samtools view -@ 32 -b -h -q 50 -N read-IDs.reads.withExc.TE-overlap.list $BAMIN | tee $BAMOUT.reads.withExc.TE-overlap.bam |
samtools view -@ 32 -q 50 -h -  | paftools.js sam2paf -L - >   $BAMOUT.reads.withExc.TE-overlap.paf ## all reads that have DEL  and overlap with a TE


# awk rearrange for preparing bed file for intersecting with EDTA intact
awk 'BEGIN{OFS="\t"; print "#" "tname","tstart","tend","qname","qstart","qend","qlen","strand","tlen","matches","lalnblock","Q","tp","mm","gn","go", "cg"}
    {print $6,$8,$9,  $1, $3,$4, $2, $5, $7,$10,$11,$12,$13,$14,$15,$16,$17}' $BAMOUT.reads.withExc.TE-overlap.paf | sortBed \
    > $BAMOUT.reads.withExc.TE-overlap.TE-format.bed # header is lost anyways # bed file for input bedtools intersect with EDTA intacts.



# Getting >500bp deletion coordinates with bamtobed
# intersecting to carry read names with it (3 additional cols compared to script_DEL reads_readcentric)
bedtools bamtobed -i $BAMOUT.reads.withExc.TE-overlap.bam -splitD  | sortBed | 
 bedtools complement -L -i stdin -g ${REF%%.fa*}*.genome |
 bedtools intersect -a - -b $BAMOUT.reads.withExc.TE-overlap.TE-format.bed -wa -wb |
 awk 'BEGIN{OFS="\t"}{if ($3 - $2 >=500 && $3 - $2 < 100000) {print $0}}' \
 > $BAMOUT.reads.withExc.TE-overlap.TE-format.largeDELcoord.bed


#intersecting with EDTA intact
cat $TEanno | # Have to remove ChrM and ptg entries first
bedtools intersect -v -a - -b $FREPEATSNS | # removing TEs that overlap with FRepeats (rDNA, slop Nstretches, organellar, telomeres) from TE intact; centromeres included though
bedtools intersect -a  $BAMOUT.reads.withExc.TE-overlap.TE-format.largeDELcoord.bed  -b - -wa -wb -loj -f 0.9 \
 >  $BAMOUT.reads.withExc.TE-overlap.TE-format.largeDELcoord.re-overlap.bed   # loj==left outer join
## Further formatting and cleaning, ; -> \t
grep -v Parent $BAMOUT.reads.withExc.TE-overlap.TE-format.largeDELcoord.re-overlap.bed |  sed 's/;[A-Za-z_]*=/\t/g'| sed 's/ID=//' |sed 's/,//g' | cut -f1-34 \
 > $BAMOUT.reads.withExc.TE-overlap.TE-format.largeDELcoord.re-overlap.clean.bed


# for  Rviz
cat $BAMOUT.reads.withExc.TE-overlap.TE-format.largeDELcoord.re-overlap.clean.bed |
 awk 'BEGIN{OFS="\t"}{print $1, $2-$5, $3-$5, "DEL", $7 "\n" $4, $5-$5,$6-$5, "READ", $7  "\n" $21, $24-$5, $25-$5, $30, $7}' | sort -k1,1 -k5,5 | uniq \
 > $BAMOUT.reads.withExc.TE-overlap.TE-format.largeDELcoord.re-overlap.clean.Rviz.bed



# getting list of mobilizing TEs
awk 'BEGIN{OFS="\t"}{ print $21, $24, $25, $30}'  $BAMOUT.reads.withExc.TE-overlap.TE-format.largeDELcoord.re-overlap.clean.bed | sort | uniq  > $BAMOUT.reads_TEs_DEL.list.bed

rm `ls *withExc.TE* | grep -v Rviz`




# conda deactivate