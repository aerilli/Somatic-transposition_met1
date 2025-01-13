#!/bin/bash

# amovilli

### TO BE RUN LIKE THIS ###
## ./Scripts/script_calling_somatic-insertions.SA.sh \
## -s  <PREFIX> \
## -b  <BAM INPUT BAM> \
## -o  <OUTDIR> \
## -t  <TE ANNOTATION FILE GFF> \
## -f  <OTHER REPEATS ANNOTATION FILE GFF> \
## -R  <REFERENCE GENOME FASTA> \


### EXAMPLE ###
## for SP in `cat ./list.SINGLE-PLANTS.txt` ; do
### ./Scripts/script_calling_somatic-insertions.SA.sh \
###  -s  ${SP} \
###  -b  ./Data/ALN/${SP}-aln2-WT-Chr_winnowmap.sortedRG.bam \
###  -o  ./Data/Somatic_events/insertions_SA \
###  -t  ./Data/Tsu-0.TE-anno.manual.gff \
###  -f  ./Data/Tsu-0.repeats.merged.gff \
###  -R  ./Reads/Tsu-0.reference.fa \
###   2>&1 | tee -a ./Data/logs/$(date '+%Y-%m-%d_%H-%M')_calling_TE_insertions_${SP}.log
## done

set -xvu
# set -o pipefail

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

# If FREPEATS are not provided
if [ -z "${FREPEATS}" ]; then
  echo -e "NotAChr\t0\t0" > ./fake-repeats.bed
  FREPEATS=./fake-repeats.bed
  FREPEATSNS=${FREPEATS}

fi


# Virtual env activation:
# containing bioawk; bedtools v2.30.0; samtools v1.16.1
## source activate /ebio/abt6_projects/met1_somatic_transpositions/conda/so_somatic-events



#### #### ( 1 ) #### ####

# Increase N-stretch annotation of 50 bp in FREPEATS, beacuse variants are often found in proximity to GAPS. 
## Centromeres are reteined because ATHILAs seem to mobilize (-> to be confirmed)
head -2 $FREPEATS \
  > $FREPEATSNS

bioawk -c gff '!/centromere/{if ( $feature == "N_stretch"){print $seqname, $source, $feature, $start-50, $end+50, $score, $filter, $strand, $group}  else {print $0}}'   $FREPEATS |
  grep Chr >> $FREPEATSNS






#### #### ( 2 ) #### ####
# extracting read IDs with supplementary aln (SA) that overlap with TEanno
## at the same time I exclude those TE annotations that overlap with the adjusted (increased NStretches) FREPEATS
grep -v -e ChrM -e ptg -e ChrC $TEanno | # Have to remove ChrM and ptg entries first
bedtools intersect -v -a - -b $FREPEATSNS | # removing TEs that overlap with <OTHER REPEATS> (rDNA, slop Nstretches, organellar, telomeres) from TE intact; centromeres included though
bedtools intersect -a $BAMIN -b - -bed |
cut -f4 | sort | uniq > $BAMOUT.IDs.TE-overlap.list
#### to increase stringency one could require a TE to overlap 99% of a read, by adding `-f 0.99 ` to `bedtools intersect -a $BAMIN -b - -bed`

samtools view -@ 32 -f 2048 -q 50 -N $BAMOUT.IDs.TE-overlap.list $BAMIN |
cut -f1 | sort | uniq  > $BAMOUT.IDs.TE-overlap.SA.list && rm $BAMOUT.IDs.TE-overlap.list

if [ ! -s $BAMOUT.IDs.TE-overlap.SA.list ]; then
  echo "probably there are no SA reads over TEs for $singleplant"
  exit 1
fi

 


#### #### ( 3 ) #### ####
# filtering reads that have SA and overlap with intact TEs using the list above
samtools view -@ 32 -b -h -q 50 -N $BAMOUT.IDs.TE-overlap.SA.list $BAMIN | bedtools intersect -v -a - -b $FREPEATSNS |
 samtools view -@ 32 -h |
 paftools.js sam2paf -L - >   $BAMOUT.TE-overlap.SA.paf ## all reads that have SA and overlap with a TE



# awk rearrange for preparing bed file for intersecting with TE annotation
awk 'BEGIN{OFS="\t"; print "#" "tname","tstart","tend","qname","qstart","qend","qlen","strand","tlen","matches","lalnblock","Q","tp","mm","gn","go", "cg"}
    {print $6,$8,$9,  $1, $3,$4, $2, $5, $7,$10,$11,$12,$13,$14,$15,$16,$17}' $BAMOUT.TE-overlap.SA.paf | sortBed \
    > $BAMOUT.TE-overlap.SA.TE-format.bed # header is lost anyways # bed file for input bedtools intersect with EDTA intacts.




#intersecting with TE annotation
grep -v -e ChrM -e ptg  $TEanno | # remove ChrM and ptg entries first
bedtools intersect -v -a - -b $FREPEATSNS | # removing TEs that overlap with <OTHER REPEATS> (rDNA, slop Nstretches, organellar, telomeres) from TE intact; centromeres included though
bedtools intersect -a $BAMOUT.TE-overlap.SA.TE-format.bed -f 0.99 -b - -wa -wb -loj > $BAMOUT.TE-overlap.SA.TE-format.re-overlap.bed   # loj==left outer join

### AND REMOVE READS THAT ARE "CONTAINED" INTO TES. AS IN THE BORDER(S) OF A READ HAS TO OVERLAP WITH READ TE BORDER
awk '{if (($21!="-1") && ($1==$18 && ($2-15<$21 && $2+15>$21)||($3-15<$22 && $3+15 >$22))){print $4}}' $BAMOUT.TE-overlap.SA.TE-format.re-overlap.bed |
 sort | uniq > $BAMOUT.TE-contained_reads.list




cat $BAMOUT.TE-overlap.SA.TE-format.re-overlap.bed |  sed 's/;[A-Za-z_]*=/\t/g'| sed 's/ID=//' |sed 's/,//g' | cut -f1-31 |
 grep -f $BAMOUT.TE-contained_reads.list |#> $BAMOUT.reads.withSA.TEintact-overlap.TE-format.re-overlap.clean.bed# awk rearrange for preparing bed file for merging
awk 'BEGIN{OFS="\t"; print "#", "RID","Rstart","Rend","RIDmerged","Rstartmerged","Rendmerged","Achr","Astart","Aend","Rlen","Astrand","Achrlen","","","Q","tp","mm","gn","go", "cg","TEchr","TEannotator","TEsuperfam","TEstart","TEend","TEscore","TEstrand","TEphase","TEID","TEname","TEclassification","TEseqont","TEidentity","TEmethod"} 
 {if ($27==""){
    print $4,$5,$6,   $1,$2,$3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26,   ".", ".", ".", ".", "."
 } else {
    print $4,$5,$6,   $1,$2,$3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31
}}' | sortBed |#    > $BAMOUT.reads.withSA.TEintact-overlap.TE-format.re-overlap.clean.merge-format.bed
bedtools merge -i stdin -d -50  -o collapse -header \
 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31  > $BAMOUT.TE-overlap.SA.TE-format.re-overlap.merged.bed

awk 'BEGIN{OFS="\t"}{if ($22 ~ /EDTA/ || $22 ~ /MANUAL/){print $0, $30, "TE"}else {print $0, "NO TE", "No TE"}}' $BAMOUT.TE-overlap.SA.TE-format.re-overlap.merged.bed  |
 tee $BAMOUT.reads.withSA.TEintact-overlap.TE-format.re-overlap.clean.merge-format.merged.Rviz.bed.tmp | cut -f1 | sort | uniq -c | awk '{if ($1==1){print $2}}' > lowquality.reads
grep -v -f lowquality.reads $BAMOUT.reads.withSA.TEintact-overlap.TE-format.re-overlap.clean.merge-format.merged.Rviz.bed.tmp > $BAMOUT.TE-overlap.SA.TE-format.re-overlap.merged.final.bed # && rm $BAMOUT.reads.withSA.TEintact-overlap.TE-format.re-overlap.clean.merge-format.merged.Rviz.bed.tmp





# getting list of mobilizing TEs
cat $BAMOUT.TE-overlap.SA.TE-format.re-overlap.merged.bed |
 grep -v -f lowquality.reads |
 awk 'BEGIN{OFS="\t"}{if ($30 ~ /^Chr/){print $30, $30} else {print $21,$24,$25,$30}}' |
  grep Chr | sed 's/\,/\n/g' |  sed 's/:/\t/' | sed 's/\.\./\t/' |  awk 'BEGIN {OFS="\t"}{if ($4==""){print $0,$1":"$2".."$3} else {print $0}}' | sort | uniq -c \
 > $BAMOUT.mobilizing-TEs.final.bed


#removing intermediates
rm `ls $BAMOUT* | grep -v final` && rm lowquality.reads
rm ./fake-repeats.bed

./Scripts/script_calling_somatic-insertions.retriever.sh $BAMOUT.TE-overlap.SA.TE-format.re-overlap.merged.final.bed 


# Virtual env deactivation
## conda deactivate